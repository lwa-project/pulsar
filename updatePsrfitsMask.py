#!/usr/bin/env python3

import os
import sys
import numpy
import argparse
from astropy.io import fits as astrofits

import lsl.common.progress as progress
from lsl.statistics import robust, kurtosis
from lsl.misc import parser as aph


def main(args):
    # Parse the command line
    if args.frequencies is not None:
        values = args.frequencies.split(',')
        
        args.frequencies = []
        for v in values:
            if v.find('-') == -1:
                args.frequencies.append( float(v) )
            else:
                v1, v2 = [float(vs) for vs in v.split('-', 1)]
                v = v1
                while v <= v2:
                    args.frequencies.append( v )
                    v += 0.1
                args.frequencies.append( v2 )
    else:
        args.frequencies = []
        
    for filename in args.filename:
        print("Working on '%s'..." % os.path.basename(filename))
        
        # Open the PRSFITS file
        hdulist = astrofits.open(filename, mode='update', memmap=True)
        
        # Figure out the integration time per sub-integration so we know how 
        # many sections to work with at a time
        nPol = hdulist[1].header['NPOL']
        nSubs = hdulist[1].header['NSBLK']
        tInt = hdulist[1].data[0][0]
        nSubsChunk = int( numpy.ceil( args.duration/tInt ) )
        print("  Polarizations: %i" % nPol)
        print("  Sub-integration time: %.3f ms" % (tInt/nSubs*1000.0,))
        print("  Sub-integrations per block: %i" % nSubs)
        print("  Block integration time: %.3f ms" % (tInt*1000.0,))
        print("  Working in chunks of %i blocks (%.3f s)" % (nSubsChunk, nSubsChunk*tInt))
        
        # Figure out the SK parameters to use
        srate = hdulist[0].header['OBSBW']*1e6
        LFFT = hdulist[1].data[0][12].size
        skM = nSubsChunk*nSubs
        skN = srate // LFFT * (tInt / nSubs)
        if nPol == 1:
            skN *= 2
        skLimits = kurtosis.get_limits(args.sk_sigma, skM, N=1.0*skN)
        print("  (p)SK M: %i" % (nSubsChunk*nSubs,))
        print("  (p)SK N: %i" % skN)
        print("  (p)SK Limits: %.4f <= valid <= %.4f" % skLimits)
        
        # Figure out what to mask for the specified frequencies and report
        toMask = []
        freq = hdulist[1].data[0][12]
        for f in args.frequencies:
            metric = numpy.abs( freq - f )
            toMaskCurrent = numpy.where( metric <= 0.05 )[0]
            toMask.extend( list(toMaskCurrent) )
        if len(toMask) > 0:
            toMask = list(set(toMask))
            toMask.sort()
            print("  Masking Channels:")
            for c in toMask:
                print("    %i -> %.3f MHz" % (c, freq[c]))
                
        # Setup the progress bar
        try:
            pbar = progress.ProgressBarPlus(max=len(hdulist[1].data)/nSubsChunk, span=58)
        except AttributeError:
            pbar = progress.ProgressBar(max=len(hdulist[1].data)/nSubsChunk, span=58)
            
        # Go!
        flagged = 0
        processed = 0
        sk = numpy.zeros((nPol, LFFT)) - 99.99
        for i in range(0, (len(hdulist[1].data)//nSubsChunk)*nSubsChunk, nSubsChunk):
            ## Load in the current block of data
            blockData = []
            blockMask = None
            for j in range(i, i+nSubsChunk):
                ### Access the correct subintegration
                subint = hdulist[1].data[j]
                
                ### Pull out various bits that we need, including:
                ###  * the weight mask
                ###  * the scale and offset values - bscl and bzero
                ###  * the actual data - data
                msk = subint[13]
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
                ### to blockData
                for k in range(nSubs):
                    d = data[:,:,k]*bscl + bzero
                    d.shape += (1,)
                    blockData.append( d )
                    
                ### Save a master mask
                try:
                    blockMask *= msk
                except TypeError:
                    blockMask = msk
                    
            blockData = numpy.concatenate(blockData, axis=2)
            
            ## Compute the S-K statistics
            for p in range(nPol):
                for l in range(LFFT):
                    sk[p,l] = kurtosis.spectral_power(blockData[p,l,:], N=1.0*skN)
                    
            ## Compute the new mask - both SK and the frequency flagging
            newMask = numpy.where( (sk < skLimits[0]) | (sk > skLimits[1]), 0.0, 1.0 )
            newMask = numpy.where( newMask.mean(axis=0) <= 0.5, 0.0, 1.0 )
            for c in toMask:
                newMask[c] *= 0.0
                
            if args.replace:
                ## Replace the existing mask
                blockMask = newMask
            else:
                ## Update the existing mask
                blockMask *= newMask
                
            ## Update file
            for j in range(i, i+nSubsChunk):
                hdulist[1].data[j][13] = blockMask
                
            ## Update the counters
            processed += LFFT
            flagged += (1.0-blockMask).sum()
            
            ## Update the progress bar and remaining time estimate
            pbar.inc()
            sys.stdout.write('  %5.1f%% %s\r' % (100.0*(1.0-blockMask).sum()/LFFT, pbar.show()))
            sys.stdout.flush()
            
        # Update the progress bar with the total time used
        sys.stdout.write('  %5.1f%% %s\n' % (100.0*flagged/processed, pbar.show()))
        sys.stdout.flush()
        
        # Done
        hdulist.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='read in a PSRFITS file and update the mask to exclude frequencies and/or time windows using a frequency mask and (pseudo) spectral kurtosis', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('filename', type=str, nargs='+', 
                        help='filename to update')
    parser.add_argument('-s', '--sk-sigma', type=aph.positive_float, default=4.0, 
                        help='(p)SK masking limit in sigma')
    parser.add_argument('-f', '--frequencies', type=str, 
                        help='comma seperated list of frequency to mask in MHz')
    parser.add_argument('-d', '--duration', type=aph.positive_float, default=10.0, 
                        help='(p)SK update interval in seconds')
    parser.add_argument('-r', '--replace', action='store_true', 
                        help='replace the current weight mask rather than augment it')
    args = parser.parse_args()
    main(args)
    
