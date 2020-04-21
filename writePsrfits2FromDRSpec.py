#!/usr/bin/env python

"""
Given a DR spectrometer file, create one of more PSRFITS file(s).
"""

# Python2 compatibility
from __future__ import print_function, division
import sys
if sys.version_info < (3,):
    range = xrange
    input = raw_input
    
import os
import sys
import numpy
import ephem
import ctypes
import argparse

import psrfits_utils.psrfits_utils as pfu

from lsl.reader.ldp import DRSpecFile
from lsl.reader import errors
import lsl.astro as astro
import lsl.common.progress as progress
from lsl.statistics import robust, kurtosis
from lsl.misc import parser as aph

from _psr import *


def resolveTarget(name):
    try:
        from urllib2 import urlopen
        from urllib import urlencode, quote_plus
    except ImportError:
        from urllib.request import urlopen
        from urllib.parse import urlencode, quote_plus
    from xml.etree import ElementTree
    
    try:
        result = urlopen('https://cdsweb.u-strasbg.fr/cgi-bin/nph-sesame/-oxp/SNV?%s' % quote_plus(name))
        tree = ElementTree.fromstring(result.read())
        target = tree.find('Target')
        service = target.find('Resolver')
        coords = service.find('jpos')
        
        serviceS = service.attrib['name'].split('=', 1)[1]
        raS, decS = coords.text.split(None, 1)
        
    except (IOError, ValueError, RuntimeError):
        raS = "---"
        decS = "---"
        serviceS = "Error"
        
    return raS, decS, serviceS


def main(args):
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
    
    # Open
    idf = DRSpecFile(args.filename)
    nFramesFile = idf.get_info('nframe')
    LFFT = idf.get_info('LFFT')
    
    # Load in basic information about the data
    srate = idf.get_info('sample_rate')
    beam = idf.get_info('beam')
    central_freq1 = idf.get_info('freq1')
    central_freq2 = idf.get_info('freq2')
    data_products = idf.get_info('data_products')
    isLinear = ('XX' in data_products) or ('YY' in data_products)
    tInt = idf.get_info('tint')
    
    # Offset, if needed
    o = 0
    if args.skip != 0.0:
       o = idf.offset(args.skip)
    nFramesFile -= int(round(o/tInt))
    
    # Sub-integration block size
    nsblk = 32
    
    ## Date
    beginDate = idf.get_info('start_time')
    beginTime = beginDate.datetime
    mjd = beginDate.mjd
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
    print("Offset: %.3f s (%.0f frames)" % (o, o/tInt))
    print("---")
    
    # Create the output PSRFITS file(s)
    pfu_out = []
    if isLinear and (not args.no_summing):
        polNames = 'I'
        nPols = 1
        def reduceEngine(x):
            y = numpy.zeros((2,x.shape[1]), dtype=numpy.float32)
            y[0,:] += x[0,:]
            y[0,:] += x[1,:]
            y[1,:] += x[2,:]
            y[1,:] += x[3,:]
            return y
    else:
        args.no_summing = True
        polNames = ''.join(data_products)
        nPols = len(data_products)
        reduceEngine = lambda x: x.astype(numpy.float32)
        
    if args.four_bit_data:
        OptimizeDataLevels = OptimizeDataLevels4Bit
    else:
        OptimizeDataLevels = OptimizeDataLevels8Bit
        
    for t in range(1, 2+1):
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
        pfo.hdr.observer = "wP2FromDRSpec.py"
        pfo.hdr.source = args.source
        pfo.hdr.fd_hand = 1
        pfo.hdr.nbits = 4 if args.four_bit_data else 8
        pfo.hdr.nsblk = nsblk
        pfo.hdr.ds_freq_fact = 1
        pfo.hdr.ds_time_fact = 1
        pfo.hdr.npol = nPols
        pfo.hdr.summed_polns = 1 if (not args.no_summing) else 0
        pfo.hdr.obs_mode = "SEARCH"
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
        
    freqBaseMHz = numpy.fft.fftshift( numpy.fft.fftfreq(LFFT, d=1.0/srate) ) / 1e6
    for i in range(len(pfu_out)):
        # Define the frequencies available in the file (in MHz)
        pfu.convert2_double_array(pfu_out[i].sub.dat_freqs, freqBaseMHz + pfu_out[i].hdr.fctr, LFFT)
        
        # Define which part of the spectra are good (1) or bad (0).  All channels
        # are good except for the two outermost.
        pfu.convert2_float_array(pfu_out[i].sub.dat_weights, numpy.ones(LFFT),  LFFT)
        pfu.set_float_value(pfu_out[i].sub.dat_weights, 0,      0)
        pfu.set_float_value(pfu_out[i].sub.dat_weights, LFFT-1, 0)
        
        # Define the data scaling (default is a scale of one and an offset of zero)
        pfu.convert2_float_array(pfu_out[i].sub.dat_offsets, numpy.zeros(LFFT*nPols), LFFT*nPols)
        pfu.convert2_float_array(pfu_out[i].sub.dat_scales,  numpy.ones(LFFT*nPols),  LFFT*nPols)
        
    # Speed things along, the data need to be processed in units of 'nsblk'.  
    # Find out how many frames that corresponds to.
    chunkSize = nsblk
    chunkTime = tInt*nsblk
    
    # Calculate the SK limites for weighting
    if (not args.no_sk_flagging) and isLinear:
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
    while True:
        ## Read in the data
        try:
            readT, t, data = idf.read(chunkTime)
            siCount += 1
        except errors.EOFError:
            break
            
        ## FFT (really promote and reshape since the data are already spectra)
        spectra = data.astype(numpy.float64)
        spectra = spectra.reshape(spectra.shape[0], -1)
        
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
        description='read in DR spectrometer files and create one or more PSRFITS file(s)', 
        epilog='NOTE:  If a source name is provided and the RA or declination is not, the script will attempt to determine these values.', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('filename', type=str, 
                        help='filename to process')
    parser.add_argument('-j', '--skip', type=aph.positive_or_zero_float, default=0.0, 
                        help='skip the specified number of seconds at the beginning of the file')
    parser.add_argument('-o', '--output', type=str, 
                        help='output file basename')
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
    