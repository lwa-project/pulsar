#!/usr/bin/env python

"""
Given several DRX files observed simultaneously with different beams, create
a collection of PSRFITS files.
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
import argparse
import traceback
from datetime import datetime

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
from lsl.misc.dedispersion import get_coherent_sample_size
from lsl.misc import parser as aph

from _psr import *


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
    args.filename.sort()
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
    
    DM = float(args.DM)
    
    startTimes = []
    nFrames = []
    for filename in args.filename:
        idf = DRXFile(filename)
        o = 0
        
        # Find out how many frame sets are in each file
        srate = idf.get_info('sample_rate')
        beampols = idf.get_info('nbeampol')
        tunepol = beampols
        nFramesFile = idf.get_info('nframe')
        
        # Offset, if needed
        o = 0
        if args.skip != 0.0:
            o = idf.offset(args.skip)
        nFramesFile -= int(o*srate/4096)*tunepol
        nFrames.append( nFramesFile // tunepol )
        
        # Get the start time of the file
        startTimes.append( idf.get_info('start_time_samples') )
        
        # Validate
        try:
            if srate != srateOld:
                raise RuntimeError("Sample rate change detected in this set of files")
        except NameError:
            srateOld = srate
            
        # Done
        idf.close()
        
    ttSkip = int(fS / srate * 4096)
    spSkip = int(fS / srate)
    frameOffsets = []
    sampleOffsets = []
    tickOffsets = []
    siCountMax = []
    for filename,startTime,nFrame in zip(args.filename, startTimes, nFrames):
        diff = max(startTimes) - startTime
        frameOffsets.append( diff // ttSkip )
        diff = diff - frameOffsets[-1]*ttSkip
        sampleOffset = diff // spSkip
        sampleOffsets.append( sampleOffset )
        if sampleOffsets[-1] == 4096:
            frameOffsets[-1] += 1
            sampleOffsets[-1] %= 4096
        if args.subsample_correction:
            tickOffsets.append( max(startTimes) - (startTime + frameOffsets[-1]*ttSkip + sampleOffsets[-1]*spSkip) )
        else:
            tickOffsets.append( 0 )
            
        nFrame = nFrame - frameOffsets[-1] - 1
        nSubints = nFrame // (nsblk * LFFT // 4096)
        siCountMax.append( nSubints )
    siCountMax = min(siCountMax)
    
    print("Proposed File Time Alignment:")
    residualOffsets = []
    for filename,startTime,frameOffset,sampleOffset,tickOffset in zip(args.filename, startTimes, frameOffsets, sampleOffsets, tickOffsets):
        tStartNow = startTime
        tStartAfter = startTime + frameOffset*ttSkip + int(sampleOffset*fS/srate) + tickOffset
        residualOffset = max(startTimes) - tStartAfter
        print("  %s with %i frames, %i samples, %i ticks" % (os.path.basename(filename), frameOffset, sampleOffset, tickOffset))
        print("    before: %i" % tStartNow)
        print("    after:  %i" % tStartAfter)
        print("      residual: %i" % residualOffset)
        
        residualOffsets.append( residualOffset )
    print("Minimum Residual: %i ticks (%.1f ns)" % (min(residualOffsets), min(residualOffsets)*(1e9/fS)))
    print("Maximum Residual: %i ticks (%.1f ns)" % (max(residualOffsets), max(residualOffsets)*(1e9/fS)))
    if not args.yes:
        out = input('=> Accept? [Y/n] ')
        if out == 'n' or out == 'N':
            sys.exit()
    else:
        print("=> Accepted via the command line")
    print(" ")
    
    # Setup the processing constraints
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
        
    for c,filename,frameOffset,sampleOffset,tickOffset in zip(range(len(args.filename)), args.filename, frameOffsets, sampleOffsets, tickOffsets):
        idf = DRXFile(filename)
            
        # Find out how many frame sets are in each file
        srate = idf.get_info('sample_rate')
        beampols = idf.get_info('nbeampol')
        tunepol = beampols
        nFramesFile = idf.get_info('nframe')
        
        # Offset, if needed
        o = 0
        if args.skip != 0.0:
            o = idf.offset(args.skip)
        nFramesFile -= int(o*srate/srate)*tunepol
        
        # Additional seek for timetag alignment across the files
        o += idf.offset(frameOffset*4096/srate)
        
        ## Date
        tStart = idf.get_info('start_time') + sampleOffset*spSkip/fS + tickOffset/fS
        beginDate = tStart.datetime
        beginTime = beginDate
        mjd = tStart.mjd
        mjd_day = int(mjd)
        mjd_sec = (mjd-mjd_day)*86400
        if args.output is None:
            args.output = "drx_%05d_%s" % (mjd_day, args.source.replace(' ', ''))
            
        ## Tuning frequencies
        central_freq1 = idf.get_info('freq1')
        central_freq2 = idf.get_info('freq2')
        beam = idf.get_info('beam')
        
        ## Coherent Dedispersion Setup
        timesPerFrame = numpy.arange(4096, dtype=numpy.float64)/srate
        spectraFreq1 = numpy.fft.fftshift( numpy.fft.fftfreq(LFFT, d=1.0/srate) ) + central_freq1
        spectraFreq2 = numpy.fft.fftshift( numpy.fft.fftfreq(LFFT, d=1.0/srate) ) + central_freq2
        
        # File summary
        print("Input Filename: %s (%i of %i)" % (filename, c+1, len(args.filename)))
        print("Date of First Frame: %s (MJD=%f)" % (str(beginDate),mjd))
        print("Tune/Pols: %i" % tunepol)
        print("Tunings: %.1f Hz, %.1f Hz" % (central_freq1, central_freq2))
        print("Sample Rate: %i Hz" % srate)
        print("Sample Time: %f s" % (LFFT/srate,))
        print("Sub-block Time: %f s" % (LFFT/srate*nsblk,))
        print("Frames: %i (%.3f s)" % (nFramesFile, 4096.0*nFramesFile / srate / tunepol))
        print("---")
        print("Using FFTW Wisdom? %s" % useWisdom)
        print("DM: %.4f pc / cm^3" % DM)
        print("Samples Needed: %i, %i to %i, %i" % (get_coherent_sample_size(central_freq1-srate/2, 1.0*srate/LFFT, DM), get_coherent_sample_size(central_freq2-srate/2, 1.0*srate/LFFT, DM), get_coherent_sample_size(central_freq1+srate/2, 1.0*srate/LFFT, DM), get_coherent_sample_size(central_freq2+srate/2, 1.0*srate/LFFT, DM)))
        
        # Parameter validation
        if get_coherent_sample_size(central_freq1-srate/2, 1.0*srate/LFFT, DM) > nsblk:
            raise RuntimeError("Too few samples for coherent dedispersion.  Considering increasing the number of channels.")
        elif get_coherent_sample_size(central_freq2-srate/2, 1.0*srate/LFFT, DM) > nsblk:
            raise RuntimeError("Too few samples for coherent dedispersion.  Considering increasing the number of channels.")
            
        # Adjust the time for the padding used for coherent dedispersion
        print("MJD shifted by %.3f ms to account for padding" %  (nsblk*LFFT/srate*1000.0,))
        beginDate = idf.get_info('start_time') + nsblk*LFFT/srate
        beginTime = beginDate.datetime
        mjd = beginDate.mjd
        
        # Create the output PSRFITS file(s)
        pfu_out = []
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
            pfo.hdr.observer = "writePsrfits2Multi.py"
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
            
            ## Coherent dedispersion information
            pfo.hdr.chan_dm = DM
            
            
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
        
        # Frequency arrays for use with the phase rotator
        freq1 = central_freq1 + numpy.fft.fftshift( numpy.fft.fftfreq(LFFT, d=1.0/srate) )
        freq2 = central_freq2 + numpy.fft.fftshift( numpy.fft.fftfreq(LFFT, d=1.0/srate) )
        
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
        pbar = progress.ProgressBarPlus(max=siCountMax, span=52)
        
        # Pre-read the first frame so that we have something to pad with, if needed
        if sampleOffset != 0:
            # Pre-read the first frame
            readT, t, dataPrev = idf.read(4096/srate)
            
        # Go!
        rdr = threading.Thread(target=reader, args=(idf, chunkTime, readerQ), kwargs={'core':0})
        rdr.setDaemon(True)
        rdr.start()
        
        # Unpack - Previous data
        incoming = getFromQueue(readerQ)
        siCount, t, rawdata = incoming
        rawSpectraPrev = PulsarEngineRaw(rawdata,LFFT)
        
        # Unpack - Current data
        incoming = getFromQueue(readerQ)
        siCount, t, rawdata = incoming
        rawSpectra = PulsarEngineRaw(rawdata, LFFT)
        
        # Main Loop
        incoming = getFromQueue(readerQ)
        while incoming[0] is not None:
            ## Unpack
            siCount, t, rawdata = incoming
            
            ## Check to see where we are
            if siCount > siCountMax:
                ### Looks like we are done, allow the reader to finish
                incoming = getFromQueue(readerQ)
                continue
                
            ## Apply the sample offset
            if sampleOffset != 0:
                try:
                    dataComb[:,:4096] = dataPrev
                except NameError:
                    dataComb = numpy.zeros((rawdata.shape[0], rawdata.shape[1]+4096), dtype=rawdata.dtype)
                    dataComb[:,:4096] = dataPrev
                dataComb[:,4096:] = rawdata
                dataPrev = dataComb[:,-4096:]
                rawdata[...] = dataComb[:,sampleOffset:sampleOffset+4096*chunkSize]
                
            ## FFT
            try:
                rawSpectraNext = PulsarEngineRaw(rawdata, LFFT, rawSpectraNext)
            except NameError:
                rawSpectraNext = PulsarEngineRaw(rawdata, LFFT)
                
            ## Apply the sub-sample offset as a phase rotation
            if tickOffset != 0:
                PhaseRotator(rawSpectra, freq1, freq2, tickOffset/fS, rawSpectra)
                
            ## S-K flagging
            flag = GenerateMask(rawSpectra)
            weight1 = numpy.where( flag[:2,:].sum(axis=0) == 0, 0, 1 ).astype(numpy.float32)
            weight2 = numpy.where( flag[2:,:].sum(axis=0) == 0, 0, 1 ).astype(numpy.float32)
            ff1 = 1.0*(LFFT - weight1.sum()) / LFFT
            ff2 = 1.0*(LFFT - weight2.sum()) / LFFT
            
            ## Dedisperse
            try:
                rawSpectraDedispersed = MultiChannelCD(rawSpectra, spectraFreq1, spectraFreq2,
                                                       1.0*srate/LFFT, DM, 
                                                       rawSpectraPrev, 
                                                       rawSpectraNext, 
                                                       rawSpectraDedispersed)
            except NameError:
                rawSpectraDedispersed = MultiChannelCD(rawSpectra, spectraFreq1, spectraFreq2,
                                                       1.0*srate/LFFT, DM, 
                                                       rawSpectraPrev, 
                                                       rawSpectraNext)
                
            ## Update the state variables used to get the CD process continuous
            rawSpectraPrev[...] = rawSpectra
            rawSpectra[...] = rawSpectraNext
            
            ## Detect power
            try:
                redData = reduceEngine(rawSpectraDedispersed, redData)
            except NameError:
                redData = reduceEngine(rawSpectraDedispersed)
                
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
        if sampleOffset != 0:
            del dataComb
        del rawSpectra
        del redData
        del bzero
        del bscale
        del bdata
        
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
        description='read in several DRX files observed simultaneously with different beams, create a collection of PSRFITS files', 
        epilog='NOTE:  If a source name is provided and the RA or declination is not, the script will attempt to determine these values.', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('DM', type=aph.positive_float, 
                        help='DM in pc cm^{-3}')
    parser.add_argument('filename', type=str, nargs='+', 
                        help='filename to process')
    parser.add_argument('-j', '--skip', type=aph.positive_or_zero_float, default=0.0, 
                        help='skip the specified number of seconds at the beginning of the file')
    parser.add_argument('-o', '--output', type=str, 
                        help='output file basename')
    parser.add_argument('-c', '--nchan', type=aph.positive_int, default=512, 
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
    parser.add_argument('-t', '--subsample-correction', action='store_true', 
                        help='enable sub-sample delay correction')
    parser.add_argument('-y', '--yes', action='store_true', 
                        help='accept the file alignment as is')
    args = parser.parse_args()
    main(args)
    
