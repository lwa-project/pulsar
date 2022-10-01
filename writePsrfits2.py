#!/usr/bin/env python

"""
Given a DRX file, create one of more PSRFITS file(s).
"""

# Python2 compatibility
from __future__ import print_function, division
try:
    input = raw_input
except NameError:
    pass
        
import os
import sys
import time
import numpy
import ctypes
import signal
import argparse
import traceback

import threading
from collections import deque
from multiprocessing import cpu_count

import psrfits_utils.psrfits_utils as pfu

from lsl.reader.ldp import DRXFile
from lsl.reader import errors
import lsl.astro as astro
import lsl.common.progress as progress
from lsl.common.dp import fS
from lsl.statistics import kurtosis
from lsl.misc import parser as aph

from _psr import *
try:
    from cupsr import *
    useCUPSR = True
except ImportError:
    useCUPSR = False


MAX_QUEUE_DEPTH = 3
readerQ = deque()


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


def reader(idf, chunkTime, outQueue, core=None, verbose=True):
    # Setup
    done = False
    siCount = 0
    
    if core is not None:
        cstatus = BindToCore(core)
        if verbose:
            print('Binding reader to core %i -> %s' % (core, cstatus))
            
    try:
        while True:
            while len(outQueue) >= MAX_QUEUE_DEPTH:
                time.sleep(0.001)
                
            ## Read in the data
            try:
                readT, t, rawdata = idf.read(chunkTime)
                siCount += 1
            except errors.EOFError:
                done = True
                break
                
            ## Add it to the queue
            outQueue.append( (siCount,t,rawdata) )
            
    except Exception as e:
        lines = traceback.format_exc()
        lines = '\x1b[2KReader Error '+lines
        print(lines,)
        
    outQueue.append( (None,done) )


def getFromQueue(queueName):
    while len(queueName) == 0:
        time.sleep(0.001)
    return queueName.popleft()


def main(args):
    # Parse command line options
    global MAX_QUEUE_DEPTH
    MAX_QUEUE_DEPTH = min([args.queue_depth, 10])
    
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
    
    # FFT length
    LFFT = args.nchan
    
    # Sub-integration block size
    nsblk = args.nsblk
    
    # Open
    idf = DRXFile(args.filename)
    
    # Load in basic information about the data
    nFramesFile = idf.get_info('nframe')
    srate = idf.get_info('sample_rate')
    beampols = idf.get_info('nbeampol')
    tunepol = beampols
    
    # Offset, if needed
    o = 0
    if args.skip != 0.0:
        o = idf.offset(args.skip)
    nFramesFile -= int(o*srate/4096)*tunepol
    
    ## Date
    beginDate = idf.get_info('start_time')
    beginTime = beginDate.datetime
    mjd = beginDate.mjd
    mjd_day = int(mjd)
    mjd_sec = (mjd-mjd_day)*86400
    if args.output is None:
        args.output = "drx_%05d_%s" % (mjd_day, args.source.replace(' ', ''))
        
    ## Tuning frequencies
    central_freq1 = idf.get_info('freq1')
    central_freq2 = idf.get_info('freq2')
    beam = idf.get_info('beam')
    
    # File summary
    print("Input Filename: %s" % args.filename)
    print("Date of First Frame: %s (MJD=%f)" % (str(beginDate),mjd))
    print("Tune/Pols: %i" % tunepol)
    print("Tunings: %.1f Hz, %.1f Hz" % (central_freq1, central_freq2))
    print("Sample Rate: %i Hz" % srate)
    print("Sample Time: %f s" % (LFFT/srate,))
    print("Sub-block Time: %f s" % (LFFT/srate*nsblk,))
    print("Frames: %i (%.3f s)" % (nFramesFile, 4096.0*nFramesFile / srate / tunepol))
    print("---")
    print("Offset: %.3f s (%i frames)" % (o, o*srate//4096*tunepol))
    print("---")
    print("Using FFTW Wisdom? %s" % useWisdom)
    print("Using cupsr instead of _psr? %s" % useCUPSR)
    
    # Create the output PSRFITS file(s)
    pfu_out = []
    if (not args.no_summing):
        polNames = 'I'
        nPols = 1
        reduceEngine = CombineToIntensity
    elif args.stokes:
        polNames = 'IQUV'
        nPols = 4
        reduceEngine = CombineToStokes
    elif args.circular:
        polNames = 'LLRR'
        nPols = 2
        reduceEngine = CombineToCircular
    else:
        polNames = 'XXYY'
        nPols = 2
        reduceEngine = CombineToLinear
        
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
        pfo.hdr.dt = LFFT / srate
        
        ## Metadata about the observation/observatory/pulsar
        pfo.hdr.observer = "writePsrfits2.py"
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
        pfo.hdr.backend = "DRX"
        pfo.hdr.project_id = "Pulsar"
        pfo.hdr.ra_str = args.ra
        pfo.hdr.dec_str = args.dec
        pfo.hdr.poln_type = "LIN" if not args.circular else "CIRC"
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
    # Find out how many frames per tuning/polarization that corresponds to.
    chunkSize = nsblk*LFFT//4096
    chunkTime = LFFT/srate*nsblk
    
    # Calculate the SK limites for weighting
    if (not args.no_sk_flagging):
        skLimits = kurtosis.get_limits(4.0, 1.0*nsblk)
        
        GenerateMask = lambda x: ComputeSKMask(x, skLimits[0], skLimits[1])
    else:
        def GenerateMask(x):
            flag = numpy.ones((4, LFFT), dtype=numpy.float32)
            flag[:,0] = 0.0
            flag[:,-1] = 0.0
            return flag
            
    # Create the progress bar so that we can keep up with the conversion.
    pbar = progress.ProgressBarPlus(max=nFramesFile//(4*chunkSize), span=52)
    
    # Go!
    rdr = threading.Thread(target=reader, args=(idf, chunkTime, readerQ), kwargs={'core':0})
    rdr.setDaemon(True)
    rdr.start()
    
    # Main Loop
    incoming = getFromQueue(readerQ)
    while incoming[0] is not None:
        ## Unpack
        siCount, t, rawdata = incoming
        
        ## FFT
        try:
            rawSpectra = PulsarEngineRaw(rawdata, LFFT, rawSpectra)
        except NameError:
            rawSpectra = PulsarEngineRaw(rawdata, LFFT)
            
        ## S-K flagging
        flag = GenerateMask(rawSpectra)
        weight1 = numpy.where( flag[:2,:].sum(axis=0) == 0, 0, 1 ).astype(numpy.float32)
        weight2 = numpy.where( flag[2:,:].sum(axis=0) == 0, 0, 1 ).astype(numpy.float32)
        ff1 = 1.0*(LFFT - weight1.sum()) / LFFT
        ff2 = 1.0*(LFFT - weight2.sum()) / LFFT
        
        ## Detect power
        try:
            redData = reduceEngine(rawSpectra, redData)
        except NameError:
            redData = reduceEngine(rawSpectra)
            
        ## Optimal data scaling
        try:
            bzero, bscale, bdata = OptimizeDataLevels(redData, LFFT, bzero, bscale, bdata)
        except NameError:
            bzero, bscale, bdata = OptimizeDataLevels(redData, LFFT)
            
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
        sys.stdout.write('%5.1f%% %5.1f%% %s %2i\r' % (ff1*100, ff2*100, pbar.show(), len(readerQ)))
        sys.stdout.flush()
        
        ## Fetch another one
        incoming = getFromQueue(readerQ)
        
    rdr.join()
    
    # Update the progress bar with the total time used but only if we have
    # reached the end of the file
    if incoming[1]:
        pbar.amount = pbar.max
    sys.stdout.write('              %s %2i\n' % (pbar.show(), len(readerQ)))
    sys.stdout.flush()
    
    # And close out the files
    for pfo in pfu_out:
        pfu.psrfits_close(pfo)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='read in DRX files and create one or more PSRFITS file(s)', 
        epilog='NOTE:  If a source name is provided and the RA or declination is not, the script will attempt to determine these values.', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('filename', type=str, 
                        help='filename to process')
    parser.add_argument('-j', '--skip', type=aph.positive_or_zero_float, default=0.0, 
                        help='skip the specified number of seconds at the beginning of the file')
    parser.add_argument('-o', '--output', type=str, 
                        help='output file basename')
    parser.add_argument('-c', '--nchan', type=aph.positive_int, default=4096, 
                        help='FFT length')
    parser.add_argument('-b', '--nsblk', type=aph.positive_int, default=4096, 
                        help='number of spetra per sub-block')
    parser.add_argument('-p', '--no-sk-flagging', action='store_true', 
                        help='disable on-the-fly SK flagging of RFI')
    parser.add_argument('-n', '--no-summing', action='store_true', 
                        help='do not sum linear polarizations')
    pgroup = parser.add_mutually_exclusive_group(required=False)
    pgroup.add_argument('-i', '--circular', action='store_true', 
                        help='convert data to RR/LL')
    pgroup.add_argument('-k', '--stokes', action='store_true', 
                        help='convert data to full Stokes')
    parser.add_argument('-s', '--source', type=str, 
                        help='source name')
    parser.add_argument('-r', '--ra', type=aph.hours, 
                        help='right ascension; HH:MM:SS.SS, J2000')
    parser.add_argument('-d', '--dec', type=aph.degrees, 
                        help='declination; sDD:MM:SS.S, J2000')
    parser.add_argument('-4', '--four-bit-data', action='store_true', 
                        help='save the spectra in 4-bit mode instead of 8-bit mode')
    parser.add_argument('-q', '--queue-depth', type=aph.positive_int, default=3, 
                        help='reader queue depth')
    args = parser.parse_args()
    main(args)
    
