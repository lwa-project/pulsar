#!/usr/bin/env python3

"""
Given a PSRFITS file, create a HDF5 file in the standard LWA1 format.
"""

import os
import re
import sys
import h5py
import numpy
import ctypes
import argparse
from datetime import datetime

from astropy.time import Time as AstroTime
from astropy.io import fits as astrofits

import data as hdfData

import lsl.astro as astro
import lsl.common.progress as progress
from lsl.misc import parser as aph


_fnRE = re.compile('.*_b(?P<beam>[1-4])(t(?P<tuning>[12]))?_.*\.fits')


def main(args):
    # Parse command line options
    filenames = args.filename
    
    # Variable that will be associated with the HDF5 instance
    f = None
    
    # Open the file and load in basic information about the observation's goal
    for c,filename in enumerate(filenames):
        ## Ready the PSRFITS file
        hdulist = astrofits.open(filename, memmap=True)
        
        ## Try to find out the beam/tuning
        mtch = _fnRE.search(filename)
        
        if mtch is None:
            beam = 0
            tuning = 1
        else:
            beam = int(mtch.group('beam'))
            try:
                tuning = int(mtch.group('tuning'))
            except:
                tuning = 1
                
        ## File specifics
        sourceName = hdulist[0].header['SRC_NAME']
        ra = hdulist[0].header['RA']
        ra = ra.split(':', 2)
        ra = sum([float(v)/60**i for i,v in enumerate(ra)])*15.0
        dec = hdulist[0].header['DEC']
        decSign = -1.0 if dec.find('-') != -1 else 1.0
        dec = dec.replace('-', '')
        dec = dec.split(':', 2)
        dec = decSign*sum([float(v)/60**i for i,v in enumerate(dec)])
        epoch = float(hdulist[0].header['EQUINOX'])
        
        tStart = AstroTime(hdulist[0].header['STT_IMJD'], (hdulist[0].header['STT_SMJD'] + hdulist[0].header['STT_OFFS'])/86400.0,
                           format='mjd', scale='utc')
        cFreq = hdulist[0].header['OBSFREQ']*1e6	# MHz -> Hz
        srate = hdulist[0].header['OBSBW']*1e6		# MHz -> Hz
        LFFT = hdulist[1].header['NCHAN']
        tInt = hdulist[1].header['TBIN']
        nSubs = hdulist[1].header['NSBLK']
        tSubs = nSubs*tInt
        nPol = hdulist[1].header['NPOL']
        if nPol == 1:
            data_products = ['I',]
        elif nPol == 2:
            if hdulist[0].header['FD_POLN'] == 'CIRC':
                data_products = ['LL', 'RR']
            else:
                data_products = ['XX', 'YY']
        else:
            data_products = ['I', 'Q', 'U', 'V']
        nChunks = len(hdulist[1].data)
        
        ## File cross-validation to make sure that everything is in order
        try:
            validationPass = True
            for keyword in ('sourceName', 'ra', 'dec', 'epoch', 'tStart', 'srate', 'LFFT', 'tInt', 'tSubs', 'nPol', 'nChunks'):
                keywordOld = keyword+"Old"
                valid = eval(f"{keyword} == {keywordOld}")
                if not valid:
                    print(f"ERROR:  Detail '{keyword}' of {os.path.basename(filename)} does not match that of the first file")
                    print("ERROR:  Aborting")
                    validationPass = False
                    
            if not validationPass:
                continue
                
        except NameError as e:
            sourceNameOld = sourceName
            raOld = ra
            decOld = dec
            epochOld = epoch
            tStartOld = tStart
            srateOld = srate
            LFFTOld = LFFT
            tIntOld = tInt
            tSubsOld = tSubs
            nPolOld = nPol
            nChunksOld = nChunks
            
        ## Pre-process the start time
        tStartI = int(tStart.unix)
        tStartF = tStart.unix - tStartI
        
        ## Convert the skip and duration values to subblocks
        skip = int(round(args.skip / tSubs))
        dur  = int(round(args.duration / tSubs))
        dur  = dur if dur else 1
        args.skip = skip * tSubs
        args.duration = dur * tSubs
        
        ## Report
        print(f"Filename: {filename} ({c+1} of {len(filenames)})")
        print(f"Date of First Frame: {tStart.datetime}")
        print(f"Beam: {beam}")
        print(f"Tuning: {tuning}")
        print(f"Sample Rate: {srate} Hz")
        print(f"Tuning Frequency: {cFreq:.3f} Hz")
        print("---")
        print(f"Target: {sourceName}")
        print(f"RA: {ra/15.0:.3f} hours")
        print(f"Dec: %.3f degrees" % dec)
        print(f"Data Products: {','.join(data_products)}")
        print(f"Integration Time: {tInt*1e3:.3f} ms")
        print(f"Sub-integrations: {nChunks} ({nChunks*tSubs:.3f} s)")
        print("---")
        print(f"Offset: {args.skip:.3f} s ({skip} subints.)")
        print(f"Duration: {args.duration:.3f} s ({dur} subints.)")
        print(f"Transform Length: {LFFT}")
        
        ## Prepare the HDF5 file
        if f is None:
            ### Create the HDF5 file
            outname = os.path.split(filename)[1]
            outname = os.path.splitext(outname)[0]
            if len(filenames) == 2:
                outname = outname.replace('t1', '')
                outname = outname.replace('t2', '')
            outname = '%s.hdf5' % outname
            
            if os.path.exists(outname):
                yn = input(f"WARNING: '{outname}' exists, overwrite? [Y/n] " % outname)
                if yn not in ('n', 'N'):
                    os.unlink(outname)
                else:
                    raise RuntimeError(f"Output file '{outname}' already exists" % outname)
                    
            ### Populate the groups
            f = hdfData.create_new_file(outname)
            hdfData.fill_minimum(f, 1, beam, srate)
            for t in (1, 2):
                hdfData.create_observation_set(f, 1, t, numpy.arange(LFFT, dtype=numpy.float64), dur*nSubs, data_products)
            f.attrs['FileGenerator'] = 'writeHDF5FromPsrfits.py'
            f.attrs['InputData'] = os.path.basename(filename)
            
            ds = {}
            ds['obs1'] = hdfData.get_observation_set(f, 1)
            ds['obs1-time'] = hdfData.get_time(f, 1)
            for t in (1, 2):
                ds[f"obs1-freq{t}"] = hdfData.get_data_set(f, 1, t, 'freq')
                for p in data_products:
                    ds[f"obs1-{p}{t}"] = hdfData.get_data_set(f, 1, t, p)
                    
            ### Add in mask information
            for t in (1, 2):
                tuningInfo = ds["obs1"].get(f"Tuning{t}", None)
                maskInfo = tuningInfo.create_group("Mask")
                for p in data_products:
                    maskInfo.create_dataset(p, ds[f"obs1-{p}{t}"].shape, 'bool')
                    ds[f"obs1-mask-{p}{t}"] = maskInfo.get(p, None)
                    
            ### Target metadata
            ds['obs1'].attrs['ObservationName'] = sourceName
            ds['obs1'].attrs['TargetName'] = sourceName
            ds['obs1'].attrs['RA'] = ra/15.0
            ds['obs1'].attrs['RA_Units'] = 'hours'
            ds['obs1'].attrs['Dec'] = dec
            ds['obs1'].attrs['Dec_Units'] = 'degrees'
            ds['obs1'].attrs['Epoch'] = epoch
            ds['obs1'].attrs['TrackingMode'] = hdulist[0].header['TRK_MODE']
            
            ### Observation metadata
            ds['obs1'].attrs['tInt'] = tInt
            ds['obs1'].attrs['tInt_Units'] = 's'
            ds['obs1'].attrs['LFFT'] = LFFT
            ds['obs1'].attrs['nChan'] = LFFT
            
            ### Create the progress bar so that we can keep up with the conversion.
            pbar = progress.ProgressBarPlus(max=len(filenames)*dur)
            
        ## Frequency information
        freq = hdulist[1].data[0][12]*1e6		# MHz -> Hz
        ds['obs1'].attrs['RBW'] = freq[1]-freq[0]
        ds['obs1'].attrs['RBW_Units'] = 'Hz'
        ds[f"obs1-freq{tuning}"][:] = freq
        
        ## Read in the data and apply what ever scaling is needed
        for i in range(skip, skip+dur):
            ### Access the correct subintegration
            subint = hdulist[1].data[i]
            
            ### Pull out various bits that we need, including:
            ###  * the start time of the subint. - tOff
            ###  * the weight mask, converted to binary - msk
            ###  * the scale and offset values - bscl and bzero
            ###  * the actual data - data
            tOff = subint[1] - subint[0]/2
            msk = numpy.where(subint[13] >= 0.5, False, True)
            bzero = subint[14]
            bscl = subint[15]
            bzero.shape = (LFFT,nPol)
            bscl.shape = (LFFT,nPol)
            bzero = bzero.T
            bscl = bscl.T
            data = subint[16]
            data.shape = (nSubs,LFFT,nPol)
            data = data.T
            
            ### Apply the scaling/offset to the data and save the results 
            ### to the HDF5 file
            for j in range(nSubs):
                k = (i-skip)*nSubs + j
                t = subint[1] + tInt*(j-nSubs//2)
                d = data[:,:,j]*bscl + bzero
                
                if c == 0:
                    ds['obs1-time'][k] = (tStartI, tStartF + t)
                for l,p in enumerate(data_products):
                    ds[f"obs1-{p}{tuning}"][k,:] = d[l,:]
                    ds[f"obs1-mask-{p}{tuning}"][k,:] = msk
                    
            ### Update the progress bar and remaining time estimate
            pbar.inc()      # pylint: disable=used-before-assignment
            sys.stdout.write('%s\r' % (pbar.show()))
            sys.stdout.flush()
            
        ## Close out the current PSRFITS file
        hdulist.close()
        
        sys.stdout.write(pbar.show()+'\n')
        sys.stdout.flush()
        
    # Done
    f.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='read in a PSRFITS file and create an HDF5 file in the standard LWA1 format', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('filename', type=str, nargs='+',
                        help='filename to process')
    parser.add_argument('-s', '--skip', type=aph.positive_or_zero_float, default=0.0,
                        help='time in seconds to skip into the file')
    parser.add_argument('-d', '--duration', type=aph.positive_float, default=10.0,
                        help='amount of time to save in seconds')
    parser.add_argument('-o', '--output', type=str,
                        help='output file basename')
    args = parser.parse_args()
    main(args)
    
