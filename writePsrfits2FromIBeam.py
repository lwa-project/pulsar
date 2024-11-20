#!/usr/bin/env python3

"""
Given a raw voltage beam file from OVRO-LWA, create a PSRFITS file.
"""

import os
import sys
import numpy
import ctypes
import struct
import argparse

import psrfits_utils.psrfits_utils as pfu

import lsl.astro as astro
from lsl.reader.base import FrameTimestamp
import lsl.common.progress as progress
from lsl.statistics import robust, kurtosis
from lsl.misc import parser as aph

from _psr import *


def read_frame_ibeam1(fh):
    hdr = fh.read(15)
    if len(hdr) < 15:
        return (None, None)
        
    hdr = struct.unpack('>BBBBBHQ', hdr)
    nbeam = hdr[3]
    nchan = hdr[2]
    data_size = nbeam*nchan*2*8
    
    data = fh.read(data_size)
    if len(data) < data_size:
        return (None, None)
    data = numpy.frombuffer(data, dtype=numpy.complex64)
    data = data.reshape(nbeam,nchan,2)
    return (hdr, data)


def read_ibeam1(fh, tInt):
    nFrame = int(round(tInt * (196e6 / (2*4096))))
    data = []
    for i in range(nFrame):
        _hdr, _data = read_frame_ibeam1(fh)
        try:
            hdr
        except NameError:
            hdr = _hdr
        _data = numpy.abs(_data)**2
        data.append(_data.transpose(2,0,1).copy())
    data = numpy.array(data).copy()
    return nFrame*(2*4096/196e6), hdr[-1]*2*4096, data[:,:,0,:]


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
    # Find out where the source is if needed
    if args.source is not None:
        if args.ra is None or args.dec is None:
            tempRA, tempDec, tempService = resolveTarget('PSR '+args.source)
            print(f"{args.source} resolved to {tempRA}, {tempDec} using '{tempService}'")
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
    fh = open(args.filename, 'rb')
    hdr, _ = read_frame_ibeam1(fh)
    frameSize = fh.tell()
    fh.seek(0)
    nFramesFile = os.path.getsize(args.filename) // frameSize
    LFFT = hdr[2]
    
    # Load in basic information about the data
    srate = 196e6/(2*4096)
    beam = 1
    central_freq1 = hdr[-2]*196e6/(2*4096) + hdr[2]*196e6/(2*4096)/2
    central_freq2 = 0.0
    data_products = ['XX', 'YY']
    isLinear = ('XX' in data_products) or ('YY' in data_products)
    tInt = 1/srate
    
    # Offset, if needed
    o = 0
    if args.skip != 0.0:
       o = int(round(o/tInt))
       fh.seek(o*frameSize)
       o = o*tInt
    nFramesFile -= int(round(o/tInt))
    
    # Sub-integration block size
    nsblk = args.nsblk
    
    ## Date
    beginDate = FrameTimestamp.from_dp_timetag(hdr[-1] * 2 * 4096)
    beginTime = beginDate.datetime
    mjd = beginDate.mjd
    mjd_day = int(mjd)
    mjd_sec = (mjd-mjd_day)*86400
    if args.output is None:
        args.output = f"drx_{mjd_day:05d}_{args.source.replace(' ', '')}"
        
    # File summary
    print(f"Input Filename: {args.filename}")
    print(f"Date of First Frame: {str(beginDate)} (MJD={mjd:f})")
    print(f"Beam: {beam}")
    print(f"Tuning: {central_freq1:.1f} Hz")
    print(f"Sample Rate: {srate} Hz")
    print(f"Sample Time: {tInt:f} s")
    print(f"Sub-block Time: {tInt * nsblk:f} s")
    print(f"Data Products: {','.join(data_products)}")
    print(f"Frames: {nFramesFile} ({tInt*nFramesFile:.3f} s)")
    print("---")
    print(f"Offset: {o:.3f} s ({o / tInt:.0f} frames)")
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
        
    for t in range(1, 2):
        ## Basic structure and bounds
        pfo = pfu.psrfits()
        pfo.basefilename = f"{args.output}_b{beam}t{t}"
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
        pfo.hdr.observer = "wP2FromIBeam.py"
        pfo.hdr.source = args.source
        pfo.hdr.fd_hand = 1
        pfo.hdr.nbits = 4 if args.four_bit_data else 8
        pfo.hdr.nsblk = nsblk
        pfo.hdr.ds_freq_fact = 1
        pfo.hdr.ds_time_fact = 1
        pfo.hdr.npol = nPols
        pfo.hdr.summed_polns = 1 if (not args.no_summing) else 0
        pfo.hdr.obs_mode = "SEARCH"
        pfo.hdr.telescope = "OVRO-LWA"
        pfo.hdr.frontend = "OVRO-LWA"
        pfo.hdr.backend = "RawVBeam"
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
        
    for i in range(len(pfu_out)):
        # Define the frequencies available in the file (in MHz)
        tfreqs = (hdr[-2] + numpy.arange(LFFT)) * (196e6/(2*4096))
        tfreqs /= 1e6
        pfu.convert2_double_array(pfu_out[i].sub.dat_freqs, tfreqs, LFFT)
        
        # Define which part of the spectra are good (1) or bad (0).  All channels
        # are good except for the two outermost.
        pfu.convert2_float_array(pfu_out[i].sub.dat_weights, numpy.ones(LFFT),  LFFT)
        
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
            flag = numpy.ones((2, LFFT), dtype=numpy.float32)
            return flag
            
    # Create the progress bar so that we can keep up with the conversion.
    pbar = progress.ProgressBarPlus(max=nFramesFile//chunkSize, span=55)
    
    # Go!
    done = False
    
    siCount = 0
    spectra = numpy.zeros((2*len(data_products), LFFT*chunkSize), dtype=numpy.float64)
    while True:
        ## Read in the data
        spectra *= 0.0
        try:
            readT, t, data = read_ibeam1(fh, chunkTime)
            spectra[0,:] = data[:,0,:].ravel()
            spectra[1,:] = data[:,1,:].ravel()
            siCount += 1
        except Exception:
            break
            
        ## FFT (really nothing since we already have what we need)
        
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
        bscale1 = bscale[:nPols,:].T.ravel()
        bdata1 = bdata[:nPols,:].T.ravel()
        
        ## Write the spectra to the PSRFITS files
        for j,sp,bz,bs,wt in zip(range(1), (bdata1,), (bzero1,), (bscale1,), (weight1,)):
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
        description='read in raw voltage beam file from OVRO-LWA create a PSRFITS file', 
        epilog='NOTE:  If a source name is provided and the RA or declination is not, the script will attempt to determine these values.', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('filename', type=str, 
                        help='filename to process')
    parser.add_argument('-j', '--skip', type=aph.positive_or_zero_float, default=0.0, 
                        help='skip the specified number of seconds at the beginning of the file')
    parser.add_argument('-o', '--output', type=str, 
                        help='output file basename')
    parser.add_argument('-b', '--nsblk', type=aph.positive_int, default=4096, 
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
    
