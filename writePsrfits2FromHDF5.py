#!/usr/bin/env python

"""
Given an HDF5 file from drspec2hdf.py, create one of more PSRFITS file(s).
"""

# Python2 compatibility
from __future__ import print_function, division
try:
    input = raw_input
except NameError:
    pass
    
import os
import sys
import h5py
import numpy
import ctypes
import argparse
from datetime import datetime

from astropy.time import Time as AstroTime

import psrfits_utils.psrfits_utils as pfu

import lsl.astro as astro
import lsl.common.progress as progress
from lsl.statistics import robust, kurtosis
from lsl.misc import parser as aph

from _psr import *


def resolveTarget(name):
    from astropy import units
    from astropy.coordinates import SkyCoord
    
    try:
        coords = SkyCoord.from_name(name)
        raS = coords.ra.to_string(unit=units.hourangle, sep=':')[:13]
        decS = coords.dec.to_string(unit=units.degree, sep=':')[:13]
        serviceS = "sesame"
    except:
        raS = "---"
        decS = "---"
        serviceS = "Error"
        
    return raS, decS, serviceS


def main(args):
    # Open the file and load in basic information about the observation's goal
    fh = h5py.File(args.filename, 'r')
    if len(fh.keys()) != 1 or 'Observation1' not in fh:
        raise RuntimeError('Only HDF5 waterfall files with a single observation, labeled "Observation1", are supported')
        
    try:
        station = fh.attrs['StationName']
    except KeyError:
        station = 'lwa1'
        
    obs1 = fh['Observation1']
    if args.source is None:
        try:
            ## Load from the observation
            sourceName = obs1.attrs['TargetName']
            
            ## Validate
            assert(sourceName != '')
            
            ## Save
            args.source = sourceName
            
        except Exception as e:
            print("WARNING: Could not load source name from file")
            
    if args.ra is None or args.dec is None:
        try:
            ## Load from the observation
            ra = obs1.attrs['RA']
            if obs1.attrs['RA_Units'] == 'degrees':
                ra /= 15.0
            dec = obs1.attrs['Dec']
            decSign = '-' if dec < 0 else '+'
            dec = abs(dec)
            
            ## Validate
            assert(ra >= 0)
            assert(ra < 24)
            assert(dec <= 90)
            
            ## Save
            args.ra = '%02d:%02d:%04.1f' % (int(ra), int(ra * 60) % 60, ra * 3600 % 60)
            args.dec = '%s%02d:%02d:%04.1f' % (decSign, int(dec), int(dec * 60) % 60, dec * 3600 % 60)
            
        except Exception as e:
            print("WARNING: Could not load source RA/dec. from file")
            
    # Find out where the source is if needed
    if args.source is not None:
        if args.ra is None or args.dec is None:
            tempRA, tempDec, tempService = resolveTarget('PSR '+args.source)
            print("%s resolved to %s, %s using '%s'" % (args.source, tempRA, tempDec, tempService))
            out = input('=> Accept? [Y/n] ')
            if out == 'n' or out == 'N':
                sys.exit()
            else:
                args.ra = tempRA
                args.dec = tempDec
                
    else:
        args.source = "None"
        
    if args.ra is None:
        args.ra = "00:00:00.00"
    if args.dec is None:
        args.dec = "+00:00:00.0"
    args.ra = str(args.ra)
    args.dec = str(args.dec)
    
    ## What's in the data?
    obs1tuning1 = obs1['Tuning1']
    try:
        obs1tuning2 = obs1['Tuning2']
    except KeyError:
        obs1tuning2 = None
    
    nFramesFile = obs1['time'].shape[0]
    srate = float(obs1.attrs['sampleRate'])
    beam = int(obs1.attrs['Beam'])
    LFFT = int(obs1.attrs['LFFT'])
    nchan = int(obs1.attrs['nChan'])
    chanOffset = LFFT - nchan		# Packing offset to deal with old HDF5 files that contain only LFFT-1 channels
    central_freq1 = obs1tuning1['freq'][LFFT//2-chanOffset]
    try:
        central_freq2 = obs1tuning2['freq'][LFFT//2-chanOffset]
    except TypeError:
        central_freq2 = 0.0
    data_products = list(obs1tuning1)
    data_products.sort()
    try:
        del data_products[ data_products.index('Saturation') ]
    except ValueError:
        pass
    try:
        del data_products[ data_products.index('freq') ]
    except ValueError:
        pass
    tInt = obs1.attrs['tInt']
    
    # Sub-integration block size
    nsblk = args.nsblk
    
    ## Date
    try:
        beginATime = AstroTime(obs1['time'][0]['int'], obs1['time'][0]['frac'],
                               format=obs1['time'].attrs['format'],
                               scale=obs1['time'].attrs['scale'])
    except (KeyError, ValueError):
        beginATime = AstroTime(obs1['time'][0], format='unix', scale='utc')
    beginDate = beginATime.utc.datetime
    beginTime = beginDate
    mjd = beginATime.utc.mjd
    mjd_day = int(mjd)
    mjd_sec = (mjd-mjd_day)*86400
    if args.output is None:
        args.output = "drx_%05d_%s" % (mjd_day, args.source.replace(' ', ''))
        
    # File summary
    print("Input Filename: %s" % args.filename)
    print("Date of First Frame: %s (MJD=%f)" % (str(beginDate),mjd))
    print("Beam: %i" % beam)
    print("Tunings: %.1f Hz, %.1f Hz" % (central_freq1, central_freq2))
    print("Sample Rate: %i Hz" % srate)
    print("Sample Time: %f s" % tInt)
    print("Sub-block Time: %f s" % (tInt*nsblk,))
    print("Data Products: %s" % ','.join(data_products))
    print("Frames: %i (%.3f s)" % (nFramesFile, tInt*nFramesFile))
    print("---")
    
    # Create the output PSRFITS file(s)
    pfu_out = []
    if 'XX' in data_products and 'YY' in data_products and (not args.no_summing):
        polNames = 'I'
        nPols = 1
        def reduceEngine(x):
            y = numpy.zeros((2,x.shape[1]), dtype=numpy.float64)
            y[0,:] += x[0,:]
            y[0,:] += x[1,:]
            y[1,:] += x[2,:]
            y[1,:] += x[3,:]
            return y.astype(numpy.float32)
    elif 'I' in data_products:
        args.no_summing = False
        polNames = 'I'
        nPols = 1
        def reduceEngine(x):
            y = numpy.zeros((2,x.shape[1]), dtype=numpy.float32)
            y[0,:] = x[0,:]
            y[1,:] = x[x.shape[0]//2,:]
            return y
    else:
        args.no_summing = True
        allowed_indices = []
        allowed_products = []
        for p,pol in enumerate(data_products):
            if pol in ('XX', 'YY'):
                allowed_indices.append(p)
                allowed_products.append(pol)
        polNames = ''.join(allowed_products)
        iPols = len(data_products)
        nPols = len(allowed_products)
        if nPols == 0:
            raise RuntimeError('No valid polarization products found: %s' % (','.join(data_products),))
        def reduceEngine(x, iPols=iPols, nPols=nPols, indicies=allowed_indices):
            y = numpy.zeros((len(allowed_products),x.shape[1]), dtype=numpy.float32)
            for i,j in enumerate(indicies):
                y[i,:] = x[j,:]
                y[nPols+i,:] = x[iPols+j]
            return y
            
    if args.four_bit_data:
        OptimizeDataLevels = OptimizeDataLevels4Bit
    else:
        OptimizeDataLevels = OptimizeDataLevels8Bit
        
    for t in range(1, 2+1):
        if t == 2 and obs1tuning2 is None:
            continue
            
        ## Basic structure and bounds
        pfo = pfu.psrfits()
        pfo.basefilename = "%s_b%it%i" % (args.output, beam, t)
        pfo.filenum = 0
        pfo.tot_rows = pfo.N = pfo.T = pfo.status = pfo.multifile = 0
        pfo.rows_per_file = 32768
        
        ## Frequency, bandwidth, and channels
        if t == 1:
            pfo.hdr.fctr=central_freq1/1e6
        else:
            pfo.hdr.fctr=central_freq2/1e6
        pfo.hdr.BW = srate/1e6
        pfo.hdr.nchan = LFFT
        pfo.hdr.df = srate/1e6/LFFT
        pfo.hdr.dt = tInt
        
        ## Metadata about the observation/observatory/pulsar
        pfo.hdr.observer = "wP2FromHDF5.py"
        pfo.hdr.source = args.source
        pfo.hdr.fd_hand = 1
        pfo.hdr.nbits = 4 if args.four_bit_data else 8
        pfo.hdr.nsblk = nsblk
        pfo.hdr.ds_freq_fact = 1
        pfo.hdr.ds_time_fact = 1
        pfo.hdr.npol = nPols
        pfo.hdr.summed_polns = 1 if (not args.no_summing) else 0
        pfo.hdr.obs_mode = "SEARCH"
        if station in ('ovro-lwa', 'ovrolwa'):
            pfo.hdr.telescope = "OVRO-LWA"
            pfo.hdr.frontend = "OVRO-LWA"
            pfo.hdr.backend = "Beamformer"
        else:
            pfo.hdr.telescope = "LWA"
            pfo.hdr.frontend = "LWA"
            pfo.hdr.backend = "DRSpectrometer"
        pfo.hdr.project_id = "Pulsar"
        pfo.hdr.ra_str = args.ra
        pfo.hdr.dec_str = args.dec
        pfo.hdr.poln_type = "LIN"
        pfo.hdr.poln_order = polNames
        pfo.hdr.date_obs = str(beginTime.strftime("%Y-%m-%dT%H:%M:%S"))     
        pfo.hdr.MJD_epoch = pfu.get_ld(mjd)
        
        ## Setup the subintegration structure
        pfo.sub.tsubint = pfo.hdr.dt*pfo.hdr.nsblk
        pfo.sub.bytes_per_subint = pfo.hdr.nchan*pfo.hdr.npol*pfo.hdr.nsblk*pfo.hdr.nbits//8
        pfo.sub.dat_freqs   = pfu.malloc_doublep(pfo.hdr.nchan*8)				# 8-bytes per double @ LFFT channels
        pfo.sub.dat_weights = pfu.malloc_floatp(pfo.hdr.nchan*4)				# 4-bytes per float @ LFFT channels
        pfo.sub.dat_offsets = pfu.malloc_floatp(pfo.hdr.nchan*pfo.hdr.npol*4)		# 4-bytes per float @ LFFT channels per pol.
        pfo.sub.dat_scales  = pfu.malloc_floatp(pfo.hdr.nchan*pfo.hdr.npol*4)		# 4-bytes per float @ LFFT channels per pol.
        if args.four_bit_data:
            pfo.sub.data = pfu.malloc_ucharp(pfo.hdr.nchan*pfo.hdr.npol*pfo.hdr.nsblk)	# 1-byte per unsigned char @ (LFFT channels x pols. x nsblk sub-integrations) samples
            pfo.sub.rawdata = pfu.malloc_ucharp(pfo.hdr.nchan*pfo.hdr.npol*pfo.hdr.nsblk//2)	# 4-bits per nibble @ (LFFT channels x pols. x nsblk sub-integrations) samples
        else:
            pfo.sub.rawdata = pfu.malloc_ucharp(pfo.hdr.nchan*pfo.hdr.npol*pfo.hdr.nsblk)	# 1-byte per unsigned char @ (LFFT channels x pols. x nsblk sub-integrations) samples
        ## Create and save it for later use
        pfu.psrfits_create(pfo)
        pfu_out.append(pfo)
        
    for i,t in enumerate((obs1tuning1, obs1tuning2)):
        if i == 1 and t is None:
            continue
            
        # Define the frequencies available in the file (in MHz) making sure to correct the array
        # if chanOffset is not zero
        tfreqs = numpy.zeros(LFFT, dtype=t['freq'].dtype)
        tfreqs[chanOffset:] = t['freq'][:]/1e6
        if chanOffset != 0:
            tfreqs[:chanOffset] = (t['freq'][0] - numpy.arange(1, chanOffset+1)[::-1]*(t['freq'][1] - t['freq'][0])) / 1e6
        pfu.convert2_double_array(pfu_out[i].sub.dat_freqs, tfreqs, LFFT)
        
        # Define which part of the spectra are good (1) or bad (0).  All channels
        # are good except for the two outermost or those that are not contained in
        # the input HDF5 file.
        pfu.convert2_float_array(pfu_out[i].sub.dat_weights, numpy.ones(LFFT),  LFFT)
        pfu.set_float_value(pfu_out[i].sub.dat_weights, 0,      0)
        pfu.set_float_value(pfu_out[i].sub.dat_weights, LFFT-1, 0)
        for j in range(chanOffset):
            pfu.set_float_value(pfu_out[i].sub.dat_weights, j, 0)
            
        # Define the data scaling (default is a scale of one and an offset of zero)
        pfu.convert2_float_array(pfu_out[i].sub.dat_offsets, numpy.zeros(LFFT*nPols), LFFT*nPols)
        pfu.convert2_float_array(pfu_out[i].sub.dat_scales,  numpy.ones(LFFT*nPols),  LFFT*nPols)
        
    # Speed things along, the data need to be processed in units of 'nsblk'.  
    # Find out how many frames that corresponds to.
    chunkSize = nsblk
    
    # Calculate the SK limites for weighting
    if (not args.no_sk_flagging) and 'XX' in data_products and 'YY' in data_products:
        skN = int(tInt*srate / LFFT)
        skLimits = kurtosis.get_limits(4.0, M=1.0*nsblk, N=1.0*skN)
        
        GenerateMask = lambda x: ComputePseudoSKMask(x, LFFT, skN, skLimits[0], skLimits[1])
    else:
        def GenerateMask(x):
            flag = numpy.ones((4, LFFT), dtype=numpy.float32)
            flag[:,0] = 0.0
            flag[:,-1] = 0.0
            return flag
            
    # Create the progress bar so that we can keep up with the conversion.
    pbar = progress.ProgressBarPlus(max=nFramesFile//chunkSize, span=55)
    
    # Go!
    done = False
    
    siCount = 0
    nSubInts = nFramesFile // chunkSize
    for i in range(nSubInts):
        ## Read in the data
        data = numpy.zeros((2*len(data_products), LFFT*chunkSize), dtype=numpy.float64)
        
        for j in range(chunkSize):
            jP = j + i*chunkSize
            nTime = obs1['time'][jP]
            try:
                nTime = nTime['int'] + nTime['frac']
            except ValueError:
                pass
                
            try:
                if nTime > oTime + 1.001*tInt:
                    # pylint: disable-next=bad-string-format-type
                    print('Warning: Time tag error in subint. %i; %.3f > %.3f + %.3f' % (siCount, nTime, oTime, tInt))
            except NameError:
                pass
            oTime = nTime
            
            k = 0
            for t in (obs1tuning1, obs1tuning2):
                if t is None:
                    continue
                    
                for p in data_products:
                    data[k, j*LFFT+chanOffset:(j+1)*LFFT] = t[p][jP,:]
                    k += 1
        siCount += 1
        
        ## Are we done yet?
        if done:
            break
            
        ## FFT
        spectra = data
        
        ## S-K flagging
        flag = GenerateMask(spectra)
        weight1 = numpy.where( flag[:2,:].sum(axis=0) == 0, 0, 1 ).astype(numpy.float32)
        weight2 = numpy.where( flag[2:,:].sum(axis=0) == 0, 0, 1 ).astype(numpy.float32)
        ff1 = 1.0*(LFFT - weight1.sum()) / LFFT
        ff2 = 1.0*(LFFT - weight2.sum()) / LFFT
        
        ## Detect power
        data = reduceEngine(spectra)
        
        ## Optimal data scaling
        bzero, bscale, bdata = OptimizeDataLevels(data, LFFT)
        
        ## Polarization mangling
        bzero1 = bzero[:nPols,:].T.ravel()
        bzero2 = bzero[nPols:,:].T.ravel()
        bscale1 = bscale[:nPols,:].T.ravel()
        bscale2 = bscale[nPols:,:].T.ravel()
        bdata1 = bdata[:nPols,:].T.ravel()
        bdata2 = bdata[nPols:,:].T.ravel()
        
        ## Write the spectra to the PSRFITS files
        for j,sp,bz,bs,wt in zip(range(2), (bdata1, bdata2), (bzero1, bzero2), (bscale1, bscale2), (weight1, weight2)):
            if j == 1 and obs1tuning2 is None:
                continue
                
            ## Time
            pfu_out[j].sub.offs = (pfu_out[j].tot_rows)*pfu_out[j].hdr.nsblk*pfu_out[j].hdr.dt+pfu_out[j].hdr.nsblk*pfu_out[j].hdr.dt/2.0
            
            ## Data
            ptr, junk = sp.__array_interface__['data']
            if args.four_bit_data:
                ctypes.memmove(int(pfu_out[j].sub.data), ptr, pfu_out[j].hdr.nchan*nPols*pfu_out[j].hdr.nsblk)
            else:
                ctypes.memmove(int(pfu_out[j].sub.rawdata), ptr, pfu_out[j].hdr.nchan*nPols*pfu_out[j].hdr.nsblk)
                
            ## Zero point
            ptr, junk = bz.__array_interface__['data']
            ctypes.memmove(int(pfu_out[j].sub.dat_offsets), ptr, pfu_out[j].hdr.nchan*nPols*4)
            
            ## Scale factor
            ptr, junk = bs.__array_interface__['data']
            ctypes.memmove(int(pfu_out[j].sub.dat_scales), ptr, pfu_out[j].hdr.nchan*nPols*4)
            
            ## SK
            ptr, junk = wt.__array_interface__['data']
            ctypes.memmove(int(pfu_out[j].sub.dat_weights), ptr, pfu_out[j].hdr.nchan*4)
            
            ## Save
            pfu.psrfits_write_subint(pfu_out[j])
            
        ## Update the progress bar and remaining time estimate
        pbar.inc()
        sys.stdout.write('%5.1f%% %5.1f%% %s\r' % (ff1*100, ff2*100, pbar.show()))
        sys.stdout.flush()
        
    # Update the progress bar with the total time used
    sys.stdout.write('              %s\n' % pbar.show())
    sys.stdout.flush()
    
    # And close out the files
    for pfo in pfu_out:
        pfu.psrfits_close(pfo)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='read in DR spectrometer/HDF5 files and create one or more PSRFITS file(s)', 
        epilog='NOTE:  If a source name is provided and the RA or declination is not, the script will attempt to determine these values.', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('filename', type=str, 
                        help='filename to process')
    parser.add_argument('-j', '--skip', type=aph.positive_or_zero_float, default=0.0, 
                        help='skip the specified number of seconds at the beginning of the file')
    parser.add_argument('-o', '--output', type=str, 
                        help='output file basename')
    parser.add_argument('-b', '--nsblk', type=aph.positive_int, default=32,
                        help='number of spetra per sub-block')
    parser.add_argument('-p', '--no-sk-flagging', action='store_true', 
                        help='disable on-the-fly SK flagging of RFI')
    parser.add_argument('-n', '--no-summing', action='store_true', 
                        help='do not sum linear polarizations')
    parser.add_argument('-s', '--source', type=str, 
                        help='source name')
    parser.add_argument('-r', '--ra', type=aph.hours, 
                        help='right ascension; HH:MM:SS.SS, J2000')
    parser.add_argument('-d', '--dec', type=aph.degrees, 
                        help='declination; sDD:MM:SS.S, J2000')
    parser.add_argument('-4', '--four-bit-data', action='store_true', 
                        help='save the spectra in 4-bit mode instead of 8-bit mode')
    args = parser.parse_args()
    main(args)
    
