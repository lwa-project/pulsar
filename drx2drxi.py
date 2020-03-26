#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Given a DRX file, create two interleaved DRX (DRXI) files, one for each tuning
"""

# Python3 compatiability
from __future__ import print_function, division
import sys
if sys.version_info > (3,):
    xrange = range
    
import os
import sys
import copy
import ephem
import struct
import argparse
from datetime import datetime
from collections import deque

from lsl.reader.ldp import DRXFile
from lsl.reader import drx, errors, buffer
from lsl import astro
from lsl.common import progress
from lsl.common.dp import fS
from lsl.misc import parser as aph


class RawDRXFrame(object):
    """
    Class to help hold and work with a raw (packed) DRX frame.
    """
    
    def __init__(self, contents):
        self.contents = bytearray(contents)
        if len(self.contents) != drx.FRAME_SIZE:
            raise errors.EOFError
        if self.contents[0] != 0xDE or self.contents[1] != 0xC0 or self.contents[2] != 0xDE or self.contents[3] != 0x5c:
            raise errors.SyncError
            
    def __getitem__(self, key):
        return self.contents[key]
        
    def __setitem__(self, key, value):
        self.contents[key] = value
        
    def parse_id(self):
        id = self.contents[4]
        id = (id & 7), ((id >> 3) & 7), ((id >> 7) & 1)
        return id
        
    @property
    def timetag(self):
        timetag = 0L
        timetag |= self.contents[16] << 56
        timetag |= self.contents[17] << 48
        timetag |= self.contents[18] << 40
        timetag |= self.contents[19] << 32
        timetag |= self.contents[20] << 24
        timetag |= self.contents[21] << 16
        timetag |= self.contents[22] <<  8
        timetag |= self.contents[23]
        return timetag
        
    @property
    def tNom(self):
        tNom = (self.contents[14] << 8) | self.contents[15]
        return tNom


class RawDRXFrameBuffer(buffer.FrameBufferBase):
    """
    A sub-type of FrameBufferBase specifically for dealing with raw (packed) DRX 
    frames.  See :class:`lsl.reader.buffer.FrameBufferBase` for a description of 
    how the buffering is implemented.
    
    Keywords:
    beams
    list of beam to expect packets for
    
    tunes
    list of tunings to expect packets for
    
    pols
    list of polarizations to expect packets for
    
    nsegments
    number of ring segments to use for the buffer (default is 20)
    
    reorder
    whether or not to reorder frames returned by get() or flush() by 
    stand/polarization (default is False)
    
    The number of segements in the ring can be converted to a buffer time in 
    seconds:
    
    +----------+--------------------------------------------------+
    |          |                 DRX Filter Code                  |
    | Segments +------+------+------+------+------+-------+-------+
    |          |  1   |  2   |  3   |  4   |  5   |  6    |  7    |
    +----------+------+------+------+------+------+-------+-------+
    |    10    | 0.16 | 0.08 | 0.04 | 0.02 | 0.01 | 0.004 | 0.002 |
    +----------+------+------+------+------+------+-------+-------+
    |    20    | 0.33 | 0.16 | 0.08 | 0.04 | 0.02 | 0.008 | 0.004 |
    +----------+------+------+------+------+------+-------+-------+
    |    30    | 0.49 | 0.25 | 0.12 | 0.06 | 0.03 | 0.013 | 0.006 |
    +----------+------+------+------+------+------+-------+-------+
    |    40    | 0.66 | 0.33 | 0.16 | 0.08 | 0.03 | 0.017 | 0.008 |
    +----------+------+------+------+------+------+-------+-------+
    |    50    | 0.82 | 0.41 | 0.20 | 0.10 | 0.04 | 0.021 | 0.010 |
    +----------+------+------+------+------+------+-------+-------+
    |   100    | 1.64 | 0.82 | 0.41 | 0.20 | 0.08 | 0.042 | 0.021 |
    +----------+------+------+------+------+------+-------+-------+
    
    """
    
    def __init__(self, beams=[], tunes=[1,2], pols=[0, 1], nsegments=20, reorder=False):
        super(RawDRXFrameBuffer, self).__init__(mode='DRX', beams=beams, tunes=tunes, pols=pols, nsegments=nsegments, reorder=reorder)
        
    def get_max_frames(self):
        """
        Calculate the maximum number of frames that we expect from 
        the setup of the observations and a list of tuples that describes
        all of the possible stand/pol combination.
        """
        
        nFrames = 0
        frameList = []
        
        nFrames = len(self.beams)*len(self.tunes)*len(self.pols)
        for beam in self.beams:
            for tune in self.tunes:
                for pol in self.pols:
                    frameList.append((beam,tune,pol))
                    
        return (nFrames, frameList)
        
    def get_figure_of_merit(self, frame):
        """
        Figure of merit for sorting frames.  For DRX it is:
            <frame timetag in ticks>
        """
        
        return frame.timetag
    
    def createFill(self, key, frameParameters):
        """
        Create a 'fill' frame of zeros using an existing good
        packet as a template.
        """

        # Get a template based on the first frame for the current buffer
        fillFrame = copy.deepcopy(self.buffer[key][0])
        
        # Get out the frame parameters and fix-up the header
        beam, tune, pol = frameParameters
        fillFrame[4] = (beam & 7) | ((tune & 7) << 3) | ((pol & 1) << 7)
        
        # Zero the data for the fill packet
        fillFrame[32:] = '\x00'*4096
        
        return fillFrame


def main(args):
    # Open the file
    idf = DRXFile(args.filename)
    
    # Load in basic information about the data
    nFramesFile = idf.get_info('nFrames')
    srate = idf.get_info('sample_rate')
    ttSkip = int(round(196e6/srate))*4096
    beam = idf.get_info('beam')
    beampols = idf.get_info('beampols')
    tunepol = beampols
    
    # Offset, if needed
    args.offset = idf.offset(args.offset)
    nFramesFile -= int(args.offset*srate/4096)*tunepol
    
    ## Date
    beginDate = ephem.Date(astro.unix_to_utcjd(idf.get_info('tStart')) - astro.DJD_OFFSET)
    beginTime = beginDate.datetime()
    mjd = astro.jd_to_mjd(astro.unix_to_utcjd(idf.get_info('tStart')))
    mjd_day = int(mjd)
    mjd_sec = (mjd-mjd_day)*86400
    
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
    print("Frames: %i (%.3f s)" % (nFramesFile, 4096.0*nFramesFile / srate / tunepol))
    
    if args.count > 0:
        nCaptures = int(args.count * srate / 4096)
    else:
        nCaptures = nFramesFile/beampols
        args.count = nCaptures * 4096 / srate
    nSkip = int(args.offset * srate / 4096 )
    
    print("Seconds to Skip:  %.2f (%i captures)" % (args.offset, nSkip))
    print("Seconds to Split: %.2f (%i captures)" % (args.count, nCaptures))
    
    outname = os.path.basename(args.filename)
    outname = os.path.splitext(outname)[0]
    print("Writing %.2f s to file '%s_b%it[12].dat'" % (nCaptures*4096/srate, outname, beam))
    
    # Ready the internal interface for file access
    fh = idf.fh
        
    # Ready the output files - one for each tune/pol
    fhOut = []
    fhOut.append( open("%s_b%it1.dat" % (outname, beam), 'wb') )
    fhOut.append( open("%s_b%it2.dat" % (outname, beam), 'wb') )
    
    try:
        pb = progress.ProgressBarPlus(max=nCaptures)
    except AttributeError:
        pb = progress.ProgressBar(max=nCaptures)
        
    newFrame = bytearray([0 for i in xrange(32+4096*2)])
    
    # Setup the buffer
    buffer = RawDRXFrameBuffer(beams=[beam,], reorder=True)
    
    # Go!
    c = 0
    eofFound = False
    while c < int(nCaptures):
        if eofFound:
            break
            
        ## Load in some frames
        rFrames = deque()
        for i in xrange(tunepol):
            try:
                rFrames.append( RawDRXFrame(fh.read(drx.FRAME_SIZE)) )
                #print(rFrames[-1].id, rFrames[-1].timetag, c, i)
            except errors.EOFError:
                eofFound = True
                buffer.append(rFrames)
                break
            except errors.SyncError:
                continue
                
        buffer.append(rFrames)
        rFrames = buffer.get()
        
        ## Continue adding frames if nothing comes out.
        if rFrames is None:
            continue
            
        ## If something comes out, process it
        for tuning in (1, 2):
            ### Load
            pair0 = rFrames[2*(tuning-1)+0]
            pair1 = rFrames[2*(tuning-1)+1]
            
            ### ID manipulation
            idX = pair0[4]
            #idY = pair1[4]
            id = (0<<7) | (1<<6) | (idX&(7<<3)) | (idX&7)
            
            ### Time tag manipulation to remove the T_NOM offset
            tNomX, timetagX = pair0.tNom, pair0.timetag
            #tNomY, timetagX = pair1.tNom, pair1.timetag
            tNom = tNomX - tNomX
            timetag = timetagX - tNomX
            
            ## Check for timetag problems
            if tuning == 1:
                try:
                    ttDiff = timetag - ttLast
                    if ttDiff != ttSkip:
                        raise RuntimeError("timetag skip at %i, %i != %i (%.1f frames)" % (c, ttDiff, ttSkip, 1.0*ttDiff/ttSkip))
                except NameError:
                    pass
                ttLast = timetag
                
            ### Build the new frame
            newFrame[0:32] = pair0[0:32]
            newFrame[32:8224:2] = pair0[32:]
            newFrame[33:8224:2] = pair1[32:]
            
            ### Update the quatities that have changed
            newFrame[4] = struct.pack('>B', id)
            newFrame[14:24] = struct.pack('>HQ', tNom, timetag)
            
            ### Save
            fhOut[tuning-1].write(newFrame)
            
        c += 1
        pb.inc(amount=1)
        if c != 0 and c % 5000 == 0:
            sys.stdout.write(pb.show()+'\r')
            sys.stdout.flush()
            
    # If we've hit the end of the file and haven't read in enough frames, 
    # flush the buffer
    if eofFound or c < int(nCaptures):
        for rFrames in buffer.flush():
            if c == int(nCaptures):
                break
                
            for tuning in (1, 2):
                ### Load
                pair0 = rFrames[2*(tuning-1)+0]
                pair1 = rFrames[2*(tuning-1)+1]
                
                ### ID manipulation
                idX = pair0[4]
                #idY = pair1[4]
                id = (0<<7) | (1<<6) | (idX&(7<<3)) | (idX&7)
                
                ### Time tag manipulation to remove the T_NOM offset
                tNomX, timetagX = pair0.tNom, pair0.timetag
                #tNomY, timetagX = pair1.tNom, pair1.timetag
                tNom = tNomX - tNomX
                timetag = timetagX - tNomX
                
                ## Check for timetag problems
                if tuning == 1:
                    try:
                        ttDiff = timetag - ttLast
                        if ttDiff != ttSkip:
                            raise RuntimeError("timetag skip at %i, %i != %i (%.1f frames)" % (c, ttDiff, ttSkip, 1.0*ttDiff/ttSkip))
                    except NameError:
                        pass
                    ttLast = timetag
                    
                ### Build the new frame
                newFrame[0:32] = pair0[0:32]
                newFrame[32:8224:2] = pair0[32:]
                newFrame[33:8224:2] = pair1[32:]
                
                ### Update the quatities that have changed
                newFrame[4] = struct.pack('>B', id)
                newFrame[14:24] = struct.pack('>HQ', tNom, timetag)
                
                ### Save
                fhOut[tuning-1].write(newFrame)
                
            c += 1
            pb.inc(amount=1)
            if c != 0 and c % 5000 == 0:
                sys.stdout.write(pb.show()+'\r')
                sys.stdout.flush()
                
    # Update the progress bar with the total time used
    pb.amount = pb.max
    sys.stdout.write(pb.show()+'\n')
    sys.stdout.flush()
    for f in fhOut:
        f.close()
        
    fh.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='convert a DRX file into two polarization-interleaved DRX files, one for each tuning', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('filename', type=str, 
                        help='filename to process')
    parser.add_argument('-c', '--count', type=aph.positive_or_zero_float, default=0.0, 
                        help='number of seconds to keep')
    parser.add_argument('-o', '--offset', type=aph.positive_or_zero_float, default=0.0, 
                        help='number of seconds to skip before splitting')
    args = parser.parse_args()
    main(args)
    
