#!/usr/bin/env python3

"""
Given a collection of single pulse files from PRESTO, plot them in an interactive way.
"""

import os
import sys
import glob
import math
import time
import numpy
import shutil
import tarfile
import argparse
import tempfile
import subprocess
from datetime import datetime
from operator import sub as opsub
from multiprocessing import Pool
from scipy.special import erf
from scipy.interpolate import interp1d
from scipy.stats import scoreatpercentile as percentile, skew, kurtosis
from scipy.signal import savgol_filter as savitzky_golay
from astropy.io import fits as astrofits

from presto.infodata import infodata
from presto.residuals import read_residuals

import lsl
from lsl import astro
from lsl.misc.dedispersion import _D, delay, incoherent
from lsl.misc.mathutils import to_dB, from_dB
from lsl.misc import parser as aph

import tkinter as tk
from tkinter import ttk, messagebox, filedialog, font as tkfont

import matplotlib
matplotlib.use('TkAgg')
matplotlib.interactive(True)

from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk, FigureCanvasTkAgg
from matplotlib.colors import Normalize
from matplotlib.collections import CircleCollection
from matplotlib import cm
from matplotlib.figure import Figure

__version__ = "0.1"
__author__ = "Jayce Dowell"


def get_aspect(ax):
    """
    Function to return the aspect ratio of a figure.

    From:
      https://stackoverflow.com/questions/41597177/get-aspect-ratio-of-axes
    """

    # Total figure size
    figW, figH = ax.get_figure().get_size_inches()
    # Axis size on figure
    _, _, w, h = ax.get_position().bounds
    # Ratio of display units
    disp_ratio = (figH * h) / (figW * w)
    # Ratio of data units
    # Negative over negative because of the order of subtraction
    data_ratio = opsub(*ax.get_ylim()) / opsub(*ax.get_xlim())

    return disp_ratio / data_ratio


def telescope2tempo(tel):
    """
    Simple function that provides the functionality of the PRESTO
    misc_utils.c telescope_to_tempocode() function.  This is needed by the
    getBarycentricCorrectionFunction() function to work with TEMPO.  Returns
    a two letter TEMPO observatory code.
    """

    tempoObservatorCodes = { 'gbt': 'GB',
                             'arecibo': 'AO',
                             'vla': 'VL',
                             'parkes': 'PK',
                             'jodrell': 'JB',
                             'gb43m': 'G1',
                             'gb 140ft': 'G1',
                             'nrao20': 'G1',
                             'nancay': 'NC',
                             'effelsberg': 'EF',
                             'srt': 'SR',
                             'wsrt': 'WT',
                             'gmrt': 'GM',
                             'lofar': 'LF',
                             'lwa': 'LW',
                             'mwa': 'MW',
                             'geocenter': 'EC', }

    try:
        return tempoObservatorCodes[tel.lower()]
    except KeyError:
        print("WARNING: Unknown telescope '%s', default to geocenter" % tel)
        return tempoObservatorCodes['geocenter']


def getBarycentricCorrectionFunction(fitsname):
    """
    Given a PSRFITS file, use TEMPO to determine an interpolating function
    that maps relative barycentric time (in seconds) to relative topocentric
    time (in seconds).  This conversion is needed for plotting the pulse
    candidates onto the waterfalls derived from the PSRFITS file.
    """

    # Open the file and read in the metadata
    hdulist = astrofits.open(fitsname, mode='readonly', memmap=True)
    ## Observatory and start time
    obs = telescope2tempo(hdulist[0].header['TELESCOP'])
    mjd = hdulist[0].header['STT_IMJD'] + (hdulist[0].header['STT_SMJD'] + hdulist[0].header['STT_OFFS'])/86400.0
    ## Pointing
    epoch = float(hdulist[0].header['EQUINOX'])
    ra = hdulist[0].header['RA']
    dec = hdulist[0].header['DEC']
    ## File length
    tInt = hdulist[1].header['TBIN']
    nSubs = hdulist[1].header['NSBLK']
    tSubs = nSubs*tInt
    nBlks = len(hdulist[1].data)
    tFile = nBlks*tSubs
    ## Done with the PSRFITS file
    hdulist.close()

    # Topocentric times to compute the barycentric times for
    # NOTE:  This includes padding at the end for the fitting in TEMPO
    topoMJD = mjd + numpy.arange(0, tFile+40.01, 20, dtype=numpy.float64)/86400.0

    # Write the TEMPO file for the conversion from topocentric to barycentric times
    with open('bary.tmp', 'w') as fh:
        fh.write("""C  Header Section
HEAD
PSR                 bary
NPRNT                  2
P0                   1.0 1
P1                   0.0
CLK            UTC(NIST)
PEPOCH           %19.13f
COORD              J2000
RA                    %s
DEC                   %s
DM                   0.0
EPHEM              DE200
C  TOA Section (uses ITAO Format)
C  First 8 columns must have + or -!
TOA\n""" % (mjd, ra, dec))
        for tMJD in topoMJD:
            fh.write("topocen+ %19.13f  0.00     0.0000  0.000000  %s\n" % (tMJD, obs))

    # Run TEMPO
    status = os.system('tempo bary.tmp > barycorr.out')
    if status != 0:
        ## This didn't work, skipping
        print("WARNING: Could not run TEMPO, skipping conversion function calculation")
        bary2topo = None

    else:
        ## Read in the barycentric times
        resids = read_residuals()
        baryMJD = resids.bary_TOA

        ## Compute relative times in seconds
        topoRel = (topoMJD - topoMJD[0])*86400.0
        baryRel = (baryMJD - baryMJD[0])*86400.0

        ## Build up the conversion function
        bary2topo = interp1d(baryRel, topoRel)

    # Cleanup
    for filename in ('bary.tmp', 'bary.par', 'barycorr.out', 'resid2.tmp', 'tempo.lis', 'matrix.tmp'):
        try:
            os.unlink(filename)
        except OSError:
            pass

    # Return
    return bary2topo


def findLimits(data, usedB=True):
    """
    Tiny function to speed up the computing of the data range for the colorbar.
    Returns a two-element list of the lowest and highest values.
    """

    dMin = data.min()
    if usedB:
        dMin = to_dB(dMin)
    if not numpy.isfinite(dMin):
        dMin = 0

    dMax = data.max()
    if usedB:
        dMax = to_dB(dMax)
    if not numpy.isfinite(dMax):
        dMax = dMin + 1

    return [dMin, dMax]


class LogNorm(Normalize):
    """
    Normalize a given value to the 0-1 range on a log scale
    """

    def __call__(self, value, clip=None):
        value = numpy.ma.asarray(value)
        mask = numpy.ma.getmaskarray(value)
        value = value.filled(self.vmax+1)
        if clip:
            numpy.clip(value, self.vmin, self.vmax)

        output = (value - self.vmin) / (self.vmax - self.vmin)
        output *= 9
        output += 1
        output = numpy.log10(output)

        output = numpy.ma.array(output, mask=mask)
        if output.shape == () and not mask:
            output = int(output)  # assume python scalar
        return output


class SqrtNorm(Normalize):
    """
    Normalize a given value to the 0-1 range on a square root scale
    """

    def __call__(self, value, clip=None):
        value = numpy.ma.asarray(value)
        mask = numpy.ma.getmaskarray(value)
        value = value.filled(self.vmax+1)
        if clip:
            numpy.clip(value, self.vmin, self.vmax)

        output = (value - self.vmin) / (self.vmax - self.vmin)
        output = numpy.sqrt(output)

        output = numpy.ma.array(output, mask=mask)
        if output.shape == () and not mask:
            output = int(output)  # assume python scalar
        return output


class SqrdNorm(Normalize):
    """
    Normalize a given value to the 0-1 range on a squared scale
    """

    def __call__(self, value, clip=None):
        value = numpy.ma.asarray(value)
        mask = numpy.ma.getmaskarray(value)
        value = value.filled(self.vmax+1)
        if clip:
            numpy.clip(value, self.vmin, self.vmax)

        output = (value - self.vmin) / (self.vmax - self.vmin)
        output = output**2

        output = numpy.ma.array(output, mask=mask)
        if output.shape == () and not mask:
            output = int(output)  # assume python scalar
        return output


class AsinhNorm(Normalize):
    """
    Normalize a given value to the 0-1 range on an inverse hyperbolic sine scale
    """

    def __call__(self, value, clip=None):
        value = numpy.ma.asarray(value)
        mask = numpy.ma.getmaskarray(value)
        value = value.filled(self.vmax+1)
        if clip:
            numpy.clip(value, self.vmin, self.vmax)

        output = (value - self.vmin) / (self.vmax - self.vmin)
        output = numpy.arcsinh(output*10.0/3.0) / numpy.arcsinh(10.0/3.0)

        output = numpy.ma.array(output, mask=mask)
        if output.shape == () and not mask:
            output = int(output)  # assume python scalar
        return output


class SinhNorm(Normalize):
    """
    Normalize a given value to the 0-1 range on an hyperbolic sine scale
    """

    def __call__(self, value, clip=None):
        value = numpy.ma.asarray(value)
        mask = numpy.ma.getmaskarray(value)
        value = value.filled(self.vmax+1)
        if clip:
            numpy.clip(value, self.vmin, self.vmax)

        output = (value - self.vmin) / (self.vmax - self.vmin)
        output = numpy.sinh(output*10.0/3.0) / numpy.sinh(10.0/3.0)

        output = numpy.ma.array(output, mask=mask)
        if output.shape == () and not mask:
            output = int(output)  # assume python scalar
        return output


class HistEqNorm(Normalize):
    """
    Normalize a given value to the 0-1 range using histogram equalization
    """

    def __call__(self, value, clip=None):
        value = numpy.ma.asarray(value)
        mask = numpy.ma.getmaskarray(value)
        value = value.filled(self.vmax+1)
        if clip:
            numpy.clip(value, self.vmin, self.vmax)

        hist, bins = numpy.histogram(value, bins=256)
        hist = numpy.insert(hist, 0, 0)
        hist = hist.cumsum() / float(hist.sum())
        histeq = interp1d(bins, hist, bounds_error=False, fill_value=0.0)
        output = histeq(value)

        output = numpy.ma.array(output, mask=mask)
        if output.shape == () and not mask:
            output = int(output)  # assume python scalar
        return output


class RefreshAwareToolbar(NavigationToolbar2Tk):
    """
    Sub-class of NavigationToolbar2Tk that includes a reference to a
    callback function that is called when the plot window is changed by the
    various buttons.
    """

    def __init__(self, canvas, parent, refreshCallback=None):
        super().__init__(canvas, parent)
        self.refreshCallback = refreshCallback

    def home(self, *args):
        super().home(*args)

        if self.refreshCallback is not None:
            self.refreshCallback()

    def forward(self, *args):
        super().forward(*args)

        if self.refreshCallback is not None:
            self.refreshCallback()

    def back(self, *args):
        super().back(*args)

        if self.refreshCallback is not None:
            self.refreshCallback()

    def release_zoom(self, event):
        super().release_zoom(event)

        if self.refreshCallback is not None:
            self.refreshCallback()

    def release_pan(self, event):
        super().release_pan(event)

        if self.refreshCallback is not None:
            self.refreshCallback()


class SinglePulse_GUI(object):
    def __init__(self, frame):
        self.frame = frame
        self.fullRes = False
        self.maxPoints = 5000

        self.filenames = []
        self.fitsname = None

        self.ax1a = None
        self.ax1b = None
        self.ax1c = None
        self.ax2 = None
        self.nBins = 21
        self.cmap = cm.get_cmap('jet')
        self.norm = Normalize
        self.plotSymbol = 'o'

        self.dataThreshold = [None, None, None]
        self.sizeProperty = 1
        self.colorProperty = 4
        self.dataWindow = [None, None, None, None]

        self._histogramCache = {}

        self.oldMarkT = None
        self.oldMarkD = None
        self.pulseClick = None

        self._mouseClickCache = {'1a':[], '1b':[], '1c':[], '2':[]}
        self._keyPressCache = {'1a':[], '1b':[], '1c':[], '2':[]}

    def loadData(self, filenames, threshold=5.0, timeRange=[0,numpy.inf], dmRange=[0,numpy.inf], widthRange=[0, numpy.inf], fitsname=None):
        print("Loading %i files with a pulse S/N threshold of %.1f" % (len(filenames), threshold))
        tStart = time.time()

        # If it a tarfile?
        if os.path.splitext(filenames[0])[1] in ('.tar.gz', '.tgz'):
            print("%6.3f s - Extracting tar file" % (time.time()-tStart,))

            if len(filenames) > 1:
                raise RuntimeError("Only one tarfile can be provided")

            self.tempdir = tempfile.mkdtemp(prefix='single-pulse-')
            tf = tarfile.open(filenames[0], mode='r:*')
            for entry in tf:
                if entry.isdir():
                    continue
                filename = os.path.basename(entry.name)
                filename = os.path.join(self.tempdir, filename)
                with open(filename, 'wb') as fh:
                    fh.write(tf.extractfile(entry.name).read())

            filenames = glob.glob(os.path.join(self.tempdir, '*'))

        # Save the filenames
        self.filenames = filenames
        self.fitsname = fitsname

        # Load the data
        print("%6.3f s - Extracting pulses" % (time.time()-tStart,))
        data = []
        meta = []
        for filename in filenames:
            ## Columns are DM, sigma, time (relative), sample (like time), and
            ## downfact (width) but are downselected to everything but sample
            try:
                newData = numpy.loadtxt(filename, dtype=numpy.float32)
                if len(newData.shape) == 2:
                    newData = newData[:,[0,1,2,3,4]]
                    data.append( newData )
                elif len(newData.shape) == 1 and newData.size > 0:
                    newData = newData.reshape(1,-1)
                    newData = newData[:,[0,1,2,3,4]]
                    data.append( newData )
            except ValueError:
                pass

            ## Metadata
            metaname = os.path.splitext(filename)[0]
            metaname = "%s.inf" % metaname
            meta.append( infodata(metaname) )
        meta = meta[0]
        data = numpy.concatenate(data)

        self.meta = meta
        self.data = numpy.ma.array(data, mask=numpy.zeros(data.shape, dtype=bool))
        self.data.data[:,4] *= 1000.0*self.meta.dt	# Convert width from samples to time in ms
        print("            -> Found %i pulses" % self.data.shape[0])

        print("%6.3f s - Applying time, DM, and width cuts" % (time.time()-tStart,))
        valid = numpy.where( (self.data[:,2] >= timeRange[0] ) & (self.data[:,2] <= timeRange[1] ) & \
                        (self.data[:,0] >= dmRange[0]   ) & (self.data[:,0] <= dmRange[1]   ) & \
                        (self.data[:,4] >= widthRange[0]) & (self.data[:,4] <= widthRange[1])    )[0]
        self.data = self.data[valid,:]
        print("            -> Downselected to %i pulses" % self.data.shape[0])

        if self.data.shape[0] == 0:
            raise RuntimeError("No pulses found after apply time and DM cuts, exiting")

        self.dmMin, self.dmMax = self.data[:,0].min(), self.data[:,0].max()
        self.snrMin, self.snrMax = self.data[:,1].min(), self.data[:,1].max()
        self.tMin, self.tMax = self.data[:,2].min(), self.data[:,2].max()
        self.widthMin, self.widthMax = self.data[:,4].min(), self.data[:,4].max()
        print("            -> DM range: %.3f to %.3f pc cm^-3" % (self.dmMin, self.dmMax))
        print("            -> S/N range: %.1f to %.1f" % (self.snrMin, self.snrMax))
        print("            -> Width range: %.3f to %.3f ms" % (self.widthMin, self.widthMax))

        print("%6.3f s - Sorting pulses in time" % (time.time()-tStart,))
        order = numpy.argsort(self.data.data[:,2])
        for i in range(self.data.shape[1]):
            self.data[:,i] = self.data[order,i]

        print("%6.3f s - Setting initial thresholds and plotting ranges" % (time.time()-tStart,))
        self.dataThreshold = [threshold, self.widthMin, self.widthMax]
        tPad = (self.tMax - self.tMin) * 0.02
        dPad = (self.dmMax - self.dmMin) * 0.02
        self.dataWindow = [self.tMin-tPad, self.tMax+tPad, self.dmMin-dPad, self.dmMax+dPad]
        print("            -> Minimum S/N: %.1f" % self.dataThreshold[0])
        print("            -> Minimum width %.3f ms" % self.dataThreshold[1])
        print("            -> Maximum width %.3f ms" % self.dataThreshold[2])
        print("            -> Time window padding: %.1f s" % tPad)
        print("            -> DM window padding: %.3f pc cm^-3" % dPad)

        print("%6.3f s - Setting default colorbar ranges" % (time.time() - tStart))
        sMin, wMin, wMax = self.dataThreshold
        tLow, tHigh, dmLow, dmHigh = self.dataWindow
        valid = numpy.where( (self.data[:,2] >= tLow ) & (self.data[:,2] <= tHigh ) & \
                        (self.data[:,0] >= dmLow) & (self.data[:,0] <= dmHigh) & \
                        (self.data[:,1] >= sMin ) & (self.data[:,4] >= wMin  ) & \
                        (self.data[:,4] <= wMax) )[0]
        self.limits = [None,]*self.data.shape[1]
        for i in range(self.data.shape[1]):
            self.limits[i] = findLimits(self.data[valid,i], usedB=False)

        if self.meta.bary and self.fitsname is not None:
            print("%6.3f s - Determining barycentric to topocentic correction factors" % (time.time()-tStart))
            self.bary2topo = getBarycentricCorrectionFunction(self.fitsname)
        else:
            self.bary2topo = None

        try:
            self.disconnect()
        except:
            pass

        try:
            shutil.rmtree(self.tempdir)
        except (AttributeError, OSError):
            pass
        del self.tempdir
        print("%6.3f s - Finished preparing data" % (time.time() - tStart))

    def getClosestPulse(self, t, dm):
        """
        Return the index of the pulse closest to the provided time and DM.
        """

        # Filter things with the right S/N and pulse width
        sLow, wLow, wHigh = self.dataThreshold
        valid = numpy.where( (self.data[:,1]>=sLow) \
                             & (self.data[:,4]>=wLow) \
                             & (self.data[:,4]<=wHigh) )[0]

        # Grab the current image scale
        aspect = get_aspect(self.ax2)

        # Find the best match after taking into account the plot aspect ratio
        d = ((self.data[valid,2]-t)/aspect)**2 + (self.data[valid,0]-dm)**2
        best = valid[numpy.argmin(d)]
        print('-> click at %.1f s, %.3f pc cm^-3 closest to pulse %i at %.1f, %.3f' % (t, dm, best, self.data[best,2], self.data[best,0]))

        return valid[numpy.argmin(d)]


    def selectTimeRange(self, t0, dm0, t1, dm1):
        """
        Given a time at DM0 and another time at DM1, select everything in time between.
        """

        fLow, fHigh = self.meta.lofreq, self.meta.lofreq + self.meta.BW

        slope = -_D*(1.0/fLow**2 - 1.0/fHigh**2)

        tCutLow = t0 + (self.dmMin-dm0)*slope
        tCutHigh = t1 + (self.dmMax-dm1)*slope

        valid1 = numpy.where( (self.data[:,2] >= tCutLow ) \
                              & (self.data[:,2] <= tCutHigh) )[0]

        deltaT1 = self.data[valid1,2] + (self.data[valid1,0]-dm0)*slope - t0
        deltaT2 = self.data[valid1,2] + (self.data[valid1,0]-dm1)*slope - t1
        valid2 = numpy.where( (deltaT1 >=0) & (deltaT2 <= 0) )[0]

        return valid1[valid2]

    def selectDMRange(self, t0, dm0, t1, dm1):
        """
        Given a time at DM0 and another time at DM1, select everything in DM between.
        """

        valid1 = numpy.where( (self.data[:,0] >= dm0) \
                              & (self.data[:,0] <= dm1) )[0]

        return valid1

    def selectTimeDMRange(self, t0, dm0, t1, dm1):
        """
        Given a time at DM0 and another time at DM1, select everything in between
        the time and DM boundaries
        """

        fLow, fHigh = self.meta.lofreq, self.meta.lofreq + self.meta.BW

        slope = -_D*(1.0/fLow**2 - 1.0/fHigh**2)

        tCutLow = t0 + (self.dmMin-dm0)*slope
        tCutHigh = t1 + (self.dmMax-dm1)*slope

        valid1 = numpy.where( (self.data[:,2] >= tCutLow ) & \
                        (self.data[:,2] <= tCutHigh) &
                        (self.data[:,0] >= dm0     ) & \
                        (self.data[:,0] <= dm1     ) )[0]

        deltaT1 = self.data[valid1,2] + (self.data[valid1,0]-dm0)*slope - t0
        deltaT2 = self.data[valid1,2] + (self.data[valid1,0]-dm1)*slope - t1
        valid2 = numpy.where( (deltaT1 >=0) & (deltaT2 <= 0) )[0]

        return valid1[valid2]

    def render(self):
        # Clean the old marks
        self.oldMarkT = None
        self.oldMarkD = None

        # Clear the old figures
        self.frame.figure1a.clf()
        self.frame.figure1b.clf()
        self.frame.figure1c.clf()
        self.frame.figure2.clf()

        self.connect()

    def draw(self, recompute=False, is_callback=False):
        """
        Draw the waterfall diagram and the total power with time.
        """

        try:
            tLowNew, tHighNew = self.ax2.get_xlim()
            dmLowNew, dmHighNew = self.ax2.get_ylim()
        except Exception as e:
            tLowNew, tHighNew = self.dataWindow[0], self.dataWindow[1]
            dmLowNew, dmHighNew = self.dataWindow[2], self.dataWindow[3]

        if tLowNew != self.dataWindow[0] or tHighNew != self.dataWindow[1]:
            self.dataWindow[0] = tLowNew
            self.dataWindow[1] = tHighNew
            recompute = True
        if dmLowNew != self.dataWindow[2] or dmHighNew != self.dataWindow[3]:
            self.dataWindow[2] = dmLowNew
            self.dataWindow[3] = dmHighNew
            recompute = True

        sMin, wMin, wMax = self.dataThreshold
        tLow, tHigh, dmLow, dmHigh = self.dataWindow
        valid = numpy.where( (self.data[:,2] >= tLow ) & (self.data[:,2] <= tHigh ) \
                             & (self.data[:,0] >= dmLow) & (self.data[:,0] <= dmHigh) \
                             & (self.data[:,1] >= sMin ) & (self.data[:,4] >= wMin  ) \
                             & (self.data[:,4] <= wMax)                               )[0]
        self.limits = [None,]*self.data.shape[1]
        for i in range(self.data.shape[1]):
            self.limits[i] = findLimits(self.data[valid,i], usedB=False)

        try:
            snrHist = self._histogramCache['snr']
            dmHist = self._histogramCache['dm']

        except KeyError:
            recompute = True

        if recompute:
            snrBins = numpy.linspace(self.dataThreshold[0], self.snrMax, self.nBins+1).astype(numpy.float32)
            dmBins = numpy.linspace(dmLow, dmHigh, self.nBins+1).astype(numpy.float32)

            try:
                from _helper import FastHistogram
                snrHist = FastHistogram(self.data[valid,1], bins=snrBins)
                dmHist = FastHistogram(self.data[valid,0], bins=dmBins)
            except ImportError:
                snrHist = numpy.histogram(self.data[valid,1], bins=snrBins)
                dmHist = numpy.histogram(self.data[valid,0], bins=dmBins)

            self._histogramCache['snr'] = snrHist
            self._histogramCache['dm'] = dmHist

        flagSNR, flagDM = False, False
        if len(self._mouseClickCache['1a']) > 0 or len(self._mouseClickCache['1b']) > 0:
            if len(self._mouseClickCache['1a']) > 0 and len(self._mouseClickCache['1b']) > 0:
                flagSNR, flagDM = True, True
                snrBounds = self._mouseClickCache['1a'][0]
                dmBounds  = self._mouseClickCache['1b'][0]
                valid2 = numpy.where( (self.data[valid,1] >= snrBounds[0]) & \
                                (self.data[valid,1] <= snrBounds[1]) & \
                                (self.data[valid,0] >= dmBounds[0]) & \
                                (self.data[valid,0] <= dmBounds[1]) )[0]

            elif len(self._mouseClickCache['1a']) > 0:
                flagSNR, flagDM = True, False
                snrBounds = self._mouseClickCache['1a'][0]
                valid2 = numpy.where( (self.data[valid,1] >= snrBounds[0]) & \
                                (self.data[valid,1] <= snrBounds[1]) )[0]

            else:
                flagSNR, flagDM = False, True
                dmBounds = self._mouseClickCache['1b'][0]
                valid2 = numpy.where( (self.data[valid,0] >= dmBounds[0]) & \
                                (self.data[valid,0] <= dmBounds[1]) )[0]

            valid2 = valid[valid2]
            alpha = 0.1

        else:
            valid2 = None
            alpha = 1.0

        # Plot 2 - Waterfall
        if len(valid) > self.maxPoints and not self.fullRes:
            decim = len(valid)//self.maxPoints
            validPlot = valid[::decim]
        else:
            validPlot = valid

        if not is_callback:
            self.frame.figure2.clf()
            self.ax2 = self.frame.figure2.gca()

            m = self.ax2.scatter(self.data[validPlot,2], self.data[validPlot,0],
                                 c=self.data[validPlot,self.colorProperty],
                                 s=self.data[validPlot,self.sizeProperty]*5,
                                 cmap=self.cmap, norm=self.norm(*self.limits[self.colorProperty]),
                                 alpha=alpha,
                                 marker=self.plotSymbol, edgecolors='face')
            try:
                cm = self.frame.figure.colorbar(m, use_gridspec=True)
            except:
                if len(self.frame.figure2.get_axes()) > 1:
                    self.frame.figure2.delaxes( self.frame.figure2.get_axes()[-1] )
                cm = self.frame.figure2.colorbar(m)
            if self.colorProperty == 0:
                cm.ax.set_ylabel('DM [pc cm$^{-3}$]')
            elif self.colorProperty == 1:
                cm.ax.set_ylabel('S/N')
            elif self.colorProperty == 2:
                cm.ax.set_ylabel('Elapsed Time [s]')
            else:
                cm.ax.set_ylabel('Width [ms]')

            if valid2 is not None:
                self.ax2.scatter(self.data[valid2,2], self.data[valid2,0],
                                 c='black',
                                 s=self.data[valid2,self.sizeProperty]*5,
                                 alpha=1.0,
                                 marker=self.plotSymbol, edgecolors='black')

            self.ax2.set_xlim((tLow,tHigh))
            self.ax2.set_ylim((dmLow,dmHigh))
            self.ax2.set_xlabel('Elapsed Time [s]')
            self.ax2.set_ylabel('DM [pc cm$^{-3}$]')

            if self.oldMarkT is not None:
                if recompute:
                    self.oldMarkT = None
                    self.oldMarkD = None
                    self.makeMark(*self.pulseClick)
                else:
                    self.ax2.plot(*self.oldMarkT, linestyle='-', marker='', color='red')
                    self.ax2.plot(*self.oldMarkD, linestyle='-', marker='', color='red')

            self.frame.figure2.tight_layout()
            self.frame.canvas2.draw()

        # Plot 1(a) - SNR histogram
        self.frame.figure1a.clf()
        self.ax1a = self.frame.figure1a.gca()
        self.ax1a.bar(snrHist[1][:-1], snrHist[0], width=snrHist[1][1]-snrHist[1][0], color='blue')
        self.ax1a.set_xlim((snrHist[1][0], snrHist[1][-1]))
        self.ax1a.set_ylim((snrHist[0].min(), snrHist[0].max()))
        self.ax1a.set_xlabel('S/N')
        self.ax1a.set_ylabel('Count')

        ## Flag a bar?
        if valid2 is not None and flagSNR:
            best = numpy.argmin( numpy.abs(snrBounds[0]-snrHist[1]) )
            self.ax1a.bar(snrBounds[0], snrHist[0][best], width=snrHist[1][1]-snrHist[1][0], color='black')

        self.frame.figure1a.tight_layout()
        self.frame.canvas1a.draw()

        # Plot 1(b) - DM histogram
        self.frame.figure1b.clf()
        self.ax1b = self.frame.figure1b.gca()
        self.ax1b.bar(dmHist[1][:-1], dmHist[0], width=dmHist[1][1]-dmHist[1][0], color='green')
        self.ax1b.set_xlim(self.ax2.get_ylim())
        self.ax1b.set_ylim((dmHist[0].min(), dmHist[0].max()))
        self.ax1b.set_xlabel('DM [pc cm$^{-3}$]')
        self.ax1b.set_ylabel('Count')

        ## Flag a bar?
        if valid2 is not None and flagDM:
            best = numpy.argmin( numpy.abs(dmBounds[0]-dmHist[1]) )
            self.ax1b.bar(dmBounds[0], dmHist[0][best], width=dmHist[1][1]-dmHist[1][0], color='black')

        self.frame.figure1b.tight_layout()
        self.frame.canvas1b.draw()

        # Plot 1(c) - DM vs. SNR
        self.frame.figure1c.clf()
        self.ax1c = self.frame.figure1c.gca()
        self.ax1c.scatter(self.data[validPlot,0], self.data[validPlot,1],
                        c=self.data[validPlot,self.colorProperty],
                        cmap=self.cmap, norm=self.norm(*self.limits[self.colorProperty]),
                        marker='+')
        self.ax1c.set_xlim(self.ax2.get_ylim())
        self.ax1c.set_ylim((self.data[valid,1].min(), self.data[valid,1].max()))
        self.ax1c.set_xlabel('DM [pc cm$^{-3}$]')
        self.ax1c.set_ylabel('S/N')

        self.frame.figure1c.tight_layout()
        self.frame.canvas1c.draw()

    def makeMark(self, clickTime, clickDM):
        if self.oldMarkT is not None:
            try:
                self.ax2.lines[-1].remove()
            except AttributeError:
                try:
                    del self.ax2.lines[-1]
                except:
                    pass
        if self.oldMarkD is not None:
            try:
                self.ax2.lines[-1].remove()
            except AttributeError:
                try:
                    del self.ax2.lines[-1]
                except:
                    pass

        fLow, fHigh = self.meta.lofreq, self.meta.lofreq + self.meta.BW

        slope = -_D*(1.0/fLow**2 - 1.0/fHigh**2)

        dm = numpy.linspace(self.dmMin, self.dmMax, 2)
        t = clickTime + (dm - clickDM)*slope
        self.oldMarkT = [t, dm]
        self.ax2.plot(*self.oldMarkT, linestyle='-', marker='', color='red')

        t = numpy.linspace(self.tMin-100, self.tMax+100, 2)
        dm = t*0 + clickDM
        self.oldMarkD = [t, dm]
        self.ax2.plot(*self.oldMarkD, linestyle='-', marker='', color='red')

        self.pulseClick = (clickTime, clickDM)

        self.frame.canvas2.draw()

    def connect(self):
        """
        Connect to all the events we need
        """

        self.cidpress1a = self.frame.figure1a.canvas.mpl_connect('button_press_event', self.on_press1a)
        self.cidpress1b = self.frame.figure1b.canvas.mpl_connect('button_press_event', self.on_press1b)
        self.cidpress1c = self.frame.figure1c.canvas.mpl_connect('button_press_event', self.on_press1c)
        self.cidpress2  = self.frame.figure2.canvas.mpl_connect('button_press_event', self.on_press2)
        self.cidkey2    = self.frame.figure2.canvas.mpl_connect('key_press_event', self.on_key2)
        self.cidmotion  = self.frame.figure2.canvas.mpl_connect('motion_notify_event', self.on_motion)

    def on_press1a(self, event):
        """
        On button press we will see if the mouse is over us and store some data
        """

        if event.inaxes:
            clickX = event.xdata
            clickY = event.ydata

            if event.button == 1:
                pass

            elif event.button == 2:
                if len(self._mouseClickCache['1a']) == 1:
                    self._mouseClickCache['1a'] = []

                bounds = [0, 1]
                for i,b in enumerate(self._histogramCache['snr'][1]):
                    if clickX >= b:
                        bounds[0] = b
                        bounds[1] = b + self._histogramCache['snr'][1][1] - self._histogramCache['snr'][1][0]

                self._mouseClickCache['1a'].append( bounds )
                self.draw()

            elif event.button == 3:
                if len(self._mouseClickCache['1a']) == 1:
                    self._mouseClickCache['1a'] = []
                    self.draw()

            else:
                pass

    def on_press1b(self, event):
        """
        On button press we will see if the mouse is over us and store some data
        """

        if event.inaxes:
            clickX = event.xdata
            clickY = event.ydata

            if event.button == 1:
                pass

            elif event.button == 2:
                if len(self._mouseClickCache['1b']) == 1:
                    self._mouseClickCache['1b'] = []

                bounds = [0, 1]
                for i,b in enumerate(self._histogramCache['dm'][1][:-1]):
                    if clickX >= b:
                        bounds[0] = b
                        bounds[1] = b + self._histogramCache['dm'][1][1] - self._histogramCache['dm'][1][0]

                self._mouseClickCache['1b'].append( bounds )
                self.draw()

            elif event.button == 3:
                if len(self._mouseClickCache['1b']) == 1:
                    self._mouseClickCache['1b'] = []
                    self.draw()

            else:
                pass

    def on_press1c(self, event):
        """
        On button press we will see if the mouse is over us and store some data
        """

        if event.inaxes:
            clickX = event.xdata
            clickY = event.ydata

            d = (clickX - self.data[:,0])**2 + (clickY - self.data[:,1])**2
            best = numpy.where( d == d.min() )[0][0]

            if event.button == 1:
                pass

            elif event.button == 2:
                ## Unmask
                print("Unmasking pulse at %.3f s, %.3f pc cm-3" % (self.data[best,2], self.data[best,0]))
                self.data.mask[best,:] = False

                self.draw(recompute=True)

            elif event.button == 3:
                ## Mask
                print("Masking pulse at %.3f s, %.3f pc cm-3" % (self.data[best,2], self.data[best,0]))
                self.data.mask[best,:] = True

                self.draw(recompute=True)

            else:
                pass

    def on_press2(self, event):
        """
        On button press we will see if the mouse is over us and store some data
        """

        if event.inaxes:
            clickX = event.xdata
            clickY = event.ydata

            best = self.getClosestPulse(clickX, clickY)

            if event.button == 1:
                if self.frame.toolbar.mode not in ('pan/zoom', 'zoom rect'):
                    self.makeMark(clickX, clickY)

            elif event.button == 2:
                ## Unmask
                print("Unmasking pulse at %.3f s, %.3f pc cm-3" % (self.data[best,2], self.data[best,0]))
                self.data.mask[best,:] = False

                self.draw(recompute=True)

            elif event.button == 3:
                ## Mask
                print("Masking pulse at %.3f s, %.3f pc cm-3" % (self.data[best,2], self.data[best,0]))
                self.data.mask[best,:] = True

                self.draw(recompute=True)

            else:
                pass

    def on_key2(self, event):
        """
        On key press we will see if the mouse is over us and store some data
        """

        if event.inaxes:
            clickX = event.xdata
            clickY = event.ydata

            best = self.getClosestPulse(clickX, clickY)

            if event.key == 'h':
                ## Help
                print("Pulse Window Keys:")
                print("  p - print the information about the underlying pulse")
                print("  e - export the current pulses to a file")
                print("  s - display a DM time slice")
                print("  w - display the PSRFITS waterfall for a pulse")
                print("  u - unmask pulses in a region of time")
                print("  m - mask pulses in a region of time")
                print("  y - unmask pulses in a region of time/DM")
                print("  n - mask pulses in a region of time/DM")
                print("  t - unmask pulses in a region of DM")
                print("  b - mask pulses in a region of DM")
                print("  h - print this help message")

            elif event.key == 'p':
                ## Print
                ### Recenter first
                self.makeMark(self.data[best,2], self.data[best,0])

                print("Time: %.3f s" % self.data.data[best,2])
                print("DM: %.3f pc cm^-3" % self.data.data[best,0])
                print("S/N: %.2f" % self.data.data[best,1])
                print("Width: %.3f ms" % self.data.data[best,4])
                print("Flagged? %s" % self.data.mask[best,0])
                print("===")

            elif event.key == 'e':
                ## Write
                outname = "plotSinglePulse.export"
                print("Saving to '%s'" % outname)

                ### Select the valid data
                sMin, wMin, wMax = self.dataThreshold
                tLow, tHigh, dmLow, dmHigh = self.dataWindow
                valid = numpy.where( (self.data[:,2] >= tLow ) & (self.data[:,2] <= tHigh ) & \
                                (self.data[:,0] >= dmLow) & (self.data[:,0] <= dmHigh) & \
                                (self.data[:,1] >= sMin ) & (self.data[:,4] >= wMin  ) & \
                                (self.data[:,4] <= wMax ) & (self.data.mask[:,0] == 0) )[0]

                ### Build a .inf file to use later
                infBase = os.path.splitext(self.filenames[0])[0]
                infBase = "%s.inf" % infBase

                with open(infBase, 'r') as ih:
                    with open('plotSinglePulse.inf', 'w') as fh:
                        for line in ih:
                            if len(line) < 3:
                                continue
                            fh.write(line)
                        fh.write("    pSP: based on template '%s'\n" % os.path.basename(infBase))
                        fh.write("    pSP: actual DM is %.3f to %.3f pc cm^-3\n" % (self.data[valid,0].min(), self.data[valid,0].max()))

                ### Build the .export (.singlepulse-like) file
                with open(outname, 'w') as fh:
                    fh.write("# DM      Sigma      Time (s)     Sample    Downfact\n")
                    for v in valid:
                        entry = (self.data[v,0], self.data[v,1], self.data[v,2], self.data[v,3], self.data[v,4])
                        fh.write("%6.4f  %5.2f  %11.4f  %6i  %6i\n" % entry)

                print("-> Done writing %i entries" % len(valid))

            elif event.key == 's':
                ## Time slice window
                ### Recenter first
                self.makeMark(self.data[best,2], self.data[best,0])

                SliceDisplay(self.frame, self.data[best,2], self.data[best,0], self.data[best,4])

            elif event.key == 'w':
                ## Waterfall window
                if self.fitsname is not None:
                    ### Recenter first
                    self.makeMark(self.data[best,2], self.data[best,0])

                    print("Time: %.3f s" % self.data.data[best,2])
                    print("DM: %.3f pc cm^-3" % self.data.data[best,0])
                    print("S/N: %.2f" % self.data.data[best,1])
                    print("Width: %.3f ms" % self.data.data[best,4])
                    print("Flagged? %s" % self.data.mask[best,0])

                    WaterfallDisplay(self.frame, self.fitsname, self.data[best,2], self.data[best,0], self.data[best,4])
                else:
                    print("No PSRFITS file specified, skipping")

            elif event.key == 'u':
                ## Mask a time range
                self._keyPressCache['2'].append( ('u', clickX, clickY) )

                if len(self._keyPressCache['2']) == 2:
                    (m0,t0,d0), (m1,t1,d1) = self._keyPressCache['2']
                    if m0 != m1:
                        del self._keyPressCache['2'][0]
                    else:
                        if t1 < t0:
                            temp = t0
                            t0 = t1
                            t1 = temp
                        print("Unmasking from %.3f to %.3f s" % (t0, t1))
                        try:
                            toMask = self.selectTimeRange(t0, d0, t1, d1)

                            self.data.mask[toMask,:] = False

                            self.draw(recompute=True)

                        except Exception as e:
                            pass

                        self._keyPressCache['2'] = []
                elif len(self._keyPressCache['2']) == 1:
                    print("Move the cursor to the other side of the time region to unmask and push 'u'")

            elif event.key == 'm':
                ## Mask a time range
                self._keyPressCache['2'].append( ('m', clickX, clickY) )

                if len(self._keyPressCache['2']) == 2:
                    (m0,t0,d0), (m1,t1,d1) = self._keyPressCache['2']
                    if m0 != m1:
                        del self._keyPressCache['2'][0]
                    else:
                        if t1 < t0:
                            temp = t0
                            t0 = t1
                            t1 = temp
                        print("Masking from %.3f to %.3f s" % (t0, t1))
                        try:
                            toMask = self.selectTimeRange(t0, d0, t1, d1)

                            self.data.mask[toMask,:] = True

                            self.draw(recompute=True)

                        except Exception as e:
                            pass

                        self._keyPressCache['2'] = []

                elif len(self._keyPressCache['2']) == 1:
                    print("Move the cursor to the other side of the time region to mask and push 'm'")

            elif event.key == 'y':
                ## Mask a time range
                self._keyPressCache['2'].append( ('y', clickX, clickY) )

                if len(self._keyPressCache['2']) == 2:
                    (m0,t0,d0), (m1,t1,d1) = self._keyPressCache['2']
                    if m0 != m1:
                        del self._keyPressCache['2'][0]
                    else:
                        if t1 < t0:
                            temp = t0
                            t0 = t1
                            t1 = temp
                        print("Unmasking from %.3f s, %.3f pc cm^-3 to %.3f s, %.3f pc cm^-3" % (t0, d0, t1, d1))
                        try:
                            toMask = self.selectTimeDMRange(t0, d0, t1, d1)

                            self.data.mask[toMask,:] = False

                            self.draw(recompute=True)

                        except Exception as e:
                            pass

                        self._keyPressCache['2'] = []
                elif len(self._keyPressCache['2']) == 1:
                    print("Move the cursor to the other corner of the time/DM region to unmask and push 'y'")

            elif event.key == 'n':
                ## Mask a time range
                self._keyPressCache['2'].append( ('n', clickX, clickY) )

                if len(self._keyPressCache['2']) == 2:
                    (m0,t0,d0), (m1,t1,d1) = self._keyPressCache['2']
                    if m0 != m1:
                        del self._keyPressCache['2'][0]
                    else:
                        if t1 < t0:
                            temp = t0
                            t0 = t1
                            t1 = temp
                        if d1 < d0:
                            temp = d0
                            d0 = d1
                            d1 = temp
                        print("Masking from %.3f s, %.3f pc cm^-3 to %.3f s, %.3f pc cm^-3" % (t0, d0, t1, d1))
                        try:
                            toMask = self.selectTimeDMRange(t0, d0, t1, d1)

                            self.data.mask[toMask,:] = True

                            self.draw(recompute=True)

                        except Exception as e:
                            pass

                        self._keyPressCache['2'] = []
                elif len(self._keyPressCache['2']) == 1:
                    print("Move the cursor to the other corner of the time/DM region to mask and push 'n'")

            elif event.key == 't':
                ## Mask a time range
                self._keyPressCache['2'].append( ('t', clickX, clickY) )

                if len(self._keyPressCache['2']) == 2:
                    (m0,t0,d0), (m1,t1,d1) = self._keyPressCache['2']
                    if m0 != m1:
                        del self._keyPressCache['2'][0]
                    else:
                        if d1 < d0:
                            temp = d0
                            d0 = d1
                            d1 = temp
                        print("Unmasking from %.3f to %.3f s" % (t0, t1))
                        try:
                            toMask = self.selectDMRange(t0, d0, t1, d1)

                            self.data.mask[toMask,:] = False

                            self.draw(recompute=True)

                        except Exception as e:
                            pass

                        self._keyPressCache['2'] = []
                elif len(self._keyPressCache['2']) == 1:
                    print("Move the cursor to the other side of the DM region to unmask and push 't'")

            elif event.key == 'b':
                ## Mask a time range
                self._keyPressCache['2'].append( ('b', clickX, clickY) )

                if len(self._keyPressCache['2']) == 2:
                    (m0,t0,d0), (m1,t1,d1) = self._keyPressCache['2']
                    if m0 != m1:
                        del self._keyPressCache['2'][0]
                    else:
                        if d1 < d0:
                            temp = d0
                            d0 = d1
                            d1 = temp
                        print("Masking from %.3f to %.3f s" % (t0, t1))
                        try:
                            toMask = self.selectDMRange(t0, d0, t1, d1)

                            self.data.mask[toMask,:] = True

                            self.draw(recompute=True)

                        except Exception as e:
                            pass

                        self._keyPressCache['2'] = []

                elif len(self._keyPressCache['2']) == 1:
                    print("Move the cursor to the other side of the DM region to mask and push 'b'")

            else:
                pass

    def on_motion(self, event):
        """
        On mouse motion display the data value under the cursor
        """

        if event.inaxes:
            clickX = event.xdata
            clickY = event.ydata

            self.frame.statusbar.config(text="%.3f s, %.3f pc cm^-3" % (clickX, clickY))
        else:
            self.frame.statusbar.config(text="")

    def disconnect(self):
        """
        Disconnect all the stored connection ids
        """

        self.frame.figure1a.canvas.mpl_disconnect(self.cidpress1a)
        self.frame.figure1b.canvas.mpl_disconnect(self.cidpress1b)
        self.frame.figure1c.canvas.mpl_disconnect(self.cidpress1c)
        self.frame.figure2.canvas.mpl_disconnect(self.cidpress2)
        self.frame.figure2.canvas.mpl_disconnect(self.cidkey2)
        self.frame.figure2.canvas.mpl_disconnect(self.cidmotion)


class MainWindow(tk.Tk):
    def __init__(self):
        super().__init__()

        self.dirname = ''
        self.filenames = []
        self.data = None
        self.examineFileButton = None
        self.examineWindow = None

        self.title("Single Pulse Viewer")
        self.geometry("1000x600")

    def render(self):
        self.initUI()
        self.initEvents()

    def initUI(self):
        # Status bar at the bottom
        self.statusbar = ttk.Label(self, text='', relief=tk.SUNKEN, anchor=tk.W)
        self.statusbar.pack(fill=tk.X, side=tk.BOTTOM)

        # Menu bar
        menubar = tk.Menu(self)
        self.config(menu=menubar)

        # File Menu
        fileMenu = tk.Menu(menubar, tearoff=0)
        fileMenu.add_command(label='Open', command=self.onOpen, accelerator='Ctrl+O')
        fileMenu.add_separator()
        fileMenu.add_command(label='Exit', command=self.onExit, accelerator='Ctrl+Q')
        menubar.add_cascade(label='File', menu=fileMenu)

        # Color Menu
        colorMenu = tk.Menu(menubar, tearoff=0)

        ## Color Mapping submenu
        self.color_value_var = tk.StringVar(value='width')
        vmapMenu = tk.Menu(colorMenu, tearoff=0)
        vmapMenu.add_radiobutton(label='DM', variable=self.color_value_var, value='dm', command=self.onColorValue)
        vmapMenu.add_radiobutton(label='S/N', variable=self.color_value_var, value='snr', command=self.onColorValue)
        vmapMenu.add_radiobutton(label='Time', variable=self.color_value_var, value='time', command=self.onColorValue)
        vmapMenu.add_radiobutton(label='Width', variable=self.color_value_var, value='width', command=self.onColorValue)
        colorMenu.add_cascade(label='Color Mapping', menu=vmapMenu)
        self.vmapMenu = vmapMenu

        ## Color Map submenu
        self.color_map_var = tk.StringVar(value='jet')
        cmapMenu = tk.Menu(colorMenu, tearoff=0)
        cmapMenu.add_radiobutton(label='Paired', variable=self.color_map_var, value='Paired', command=self.onColorMap)
        cmapMenu.add_radiobutton(label='Spectral', variable=self.color_map_var, value='Spectral', command=self.onColorMap)
        cmapMenu.add_radiobutton(label='Bone', variable=self.color_map_var, value='bone', command=self.onColorMap)
        cmapMenu.add_radiobutton(label='Jet', variable=self.color_map_var, value='jet', command=self.onColorMap)
        cmapMenu.add_radiobutton(label='Earth', variable=self.color_map_var, value='gist_earth', command=self.onColorMap)
        cmapMenu.add_radiobutton(label='Heat', variable=self.color_map_var, value='gist_heat', command=self.onColorMap)
        cmapMenu.add_radiobutton(label='NCAR', variable=self.color_map_var, value='gist_ncar', command=self.onColorMap)
        cmapMenu.add_radiobutton(label='Rainbow', variable=self.color_map_var, value='gist_rainbow', command=self.onColorMap)
        cmapMenu.add_radiobutton(label='Stern', variable=self.color_map_var, value='gist_stern', command=self.onColorMap)
        cmapMenu.add_radiobutton(label='Gray', variable=self.color_map_var, value='gist_gray', command=self.onColorMap)
        cmapMenu.add_separator()
        self.color_invert_var = tk.BooleanVar(value=False)
        cmapMenu.add_checkbutton(label='Invert', variable=self.color_invert_var, command=self.onColorMap)
        colorMenu.add_cascade(label='Color Map', menu=cmapMenu)
        self.cmapMenu = cmapMenu

        ## Color Stretch submenu
        self.color_stretch_var = tk.StringVar(value='linear')
        smapMenu = tk.Menu(colorMenu, tearoff=0)
        smapMenu.add_radiobutton(label='Linear', variable=self.color_stretch_var, value='linear', command=self.onColorStretch)
        smapMenu.add_radiobutton(label='Log', variable=self.color_stretch_var, value='log', command=self.onColorStretch)
        smapMenu.add_radiobutton(label='Square Root', variable=self.color_stretch_var, value='sqrt', command=self.onColorStretch)
        smapMenu.add_radiobutton(label='Squared', variable=self.color_stretch_var, value='sqrd', command=self.onColorStretch)
        smapMenu.add_radiobutton(label='ASinh', variable=self.color_stretch_var, value='asinh', command=self.onColorStretch)
        smapMenu.add_radiobutton(label='Sinh', variable=self.color_stretch_var, value='sinh', command=self.onColorStretch)
        smapMenu.add_radiobutton(label='Histogram Equalization', variable=self.color_stretch_var, value='hist', command=self.onColorStretch)
        colorMenu.add_cascade(label='Color Stretch', menu=smapMenu)
        self.smapMenu = smapMenu

        menubar.add_cascade(label='Color', menu=colorMenu)

        # Data Menu
        dataMenu = tk.Menu(menubar, tearoff=0)

        ## Plot Symbol submenu
        self.plot_symbol_var = tk.StringVar(value='o')
        mmapMenu = tk.Menu(dataMenu, tearoff=0)
        mmapMenu.add_radiobutton(label='Circle', variable=self.plot_symbol_var, value='o', command=self.onDataSymbol)
        mmapMenu.add_radiobutton(label='Square', variable=self.plot_symbol_var, value='s', command=self.onDataSymbol)
        mmapMenu.add_radiobutton(label='Diamond', variable=self.plot_symbol_var, value='D', command=self.onDataSymbol)
        mmapMenu.add_radiobutton(label='Hexagon', variable=self.plot_symbol_var, value='h', command=self.onDataSymbol)
        mmapMenu.add_radiobutton(label='Plus Sign', variable=self.plot_symbol_var, value='+', command=self.onDataSymbol)
        dataMenu.add_cascade(label='Plot Symbol', menu=mmapMenu)
        self.mmapMenu = mmapMenu

        ## Size Mapping submenu
        self.size_map_var = tk.StringVar(value='snr')
        amapMenu = tk.Menu(dataMenu, tearoff=0)
        amapMenu.add_radiobutton(label='DM', variable=self.size_map_var, value='dm', command=self.onDataSize)
        amapMenu.add_radiobutton(label='S/N', variable=self.size_map_var, value='snr', command=self.onDataSize)
        amapMenu.add_radiobutton(label='Time', variable=self.size_map_var, value='time', command=self.onDataSize)
        amapMenu.add_radiobutton(label='Width', variable=self.size_map_var, value='width', command=self.onDataSize)
        dataMenu.add_cascade(label='Size Mapping', menu=amapMenu)
        self.amapMenu = amapMenu

        dataMenu.add_command(label='Adjust Thresholds', command=self.onDataAdjust)

        menubar.add_cascade(label='Data', menu=dataMenu)

        # Display Menu
        displayMenu = tk.Menu(menubar, tearoff=0)
        self.decimate_var = tk.BooleanVar(value=not self.data.fullRes if self.data else True)
        displayMenu.add_checkbutton(label='Decimation', variable=self.decimate_var, command=self.onDisplayDecimate)
        displayMenu.add_command(label='Decimation Adjust', command=self.onDisplayDecimateAdjust)
        menubar.add_cascade(label='Display', menu=displayMenu)
        self.displayMenu = displayMenu

        # Help Menu
        helpMenu = tk.Menu(menubar, tearoff=0)
        helpMenu.add_command(label='plotSinglePulse Handbook', command=self.onHelp, accelerator='F1')
        helpMenu.add_separator()
        helpMenu.add_command(label='About', command=self.onAbout)
        menubar.add_cascade(label='Help', menu=helpMenu)

        # Main content area
        # Top panel with three histogram plots
        panel1 = ttk.Frame(self)
        panel1.pack(fill=tk.BOTH, expand=True)

        # SNR histogram
        self.figure1a = Figure(figsize=(2, 2))
        self.canvas1a = FigureCanvasTkAgg(self.figure1a, master=panel1)
        self.canvas1a.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        # DM histogram
        self.figure1b = Figure(figsize=(2, 2))
        self.canvas1b = FigureCanvasTkAgg(self.figure1b, master=panel1)
        self.canvas1b.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        # DM vs SNR scatter plot
        self.figure1c = Figure(figsize=(2, 2))
        self.canvas1c = FigureCanvasTkAgg(self.figure1c, master=panel1)
        self.canvas1c.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        # Bottom panel with main scatter plot and toolbar
        panel3 = ttk.Frame(self)
        panel3.pack(fill=tk.BOTH, expand=True)

        self.figure2 = Figure(figsize=(6, 2))
        self.canvas2 = FigureCanvasTkAgg(self.figure2, master=panel3)
        self.canvas2.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        self.toolbar = RefreshAwareToolbar(self.canvas2, panel3,
                                           refreshCallback=lambda: self.data.draw(is_callback=True))
        self.toolbar.update()
        self.panel3 = panel3

    def initEvents(self):
        # Keyboard shortcuts
        self.bind('<Control-o>', lambda e: self.onOpen())
        self.bind('<Control-q>', lambda e: self.onExit())
        self.bind('<F1>', lambda e: self.onHelp())

        # Key events on canvases
        self.canvas1a.get_tk_widget().bind('<KeyRelease>', self.onKeyPress)
        self.canvas1b.get_tk_widget().bind('<KeyRelease>', self.onKeyPress)
        self.canvas2.get_tk_widget().bind('<KeyRelease>', self.onKeyPress)

        # Window resize
        self.bind('<Configure>', self.onSize)

        # Window close
        self.protocol("WM_DELETE_WINDOW", self.onExit)

    def onOpen(self):
        """Open a file."""
        filename = filedialog.askopenfilename(
            title="Choose a file",
            initialdir=self.dirname,
            filetypes=[('Single Pulse', '*.singlepulse'), ('All Files', '*.*')]
        )
        if filename:
            self.filename = os.path.basename(filename)
            self.dirname = os.path.dirname(filename)
            self.data = SinglePulse_GUI(self)

    def onHelp(self):
        """Display the help window."""
        HelpWindow(self)

    def onAbout(self):
        """Display a brief 'about' window."""
        messagebox.showinfo('About plotSinglePulse',
            f'plotSinglePulse v{__version__}\n\n'
            f'GUI for plotting single pulse data from PRESTO.\n\n'
            f'LSL Version: {lsl.version.version}\n\n'
            f'Developer: {__author__}')

    def onExit(self):
        """Quit the application."""
        self.destroy()

    def onColorValue(self):
        """Set the color coding to the specified quantity and refresh the plots."""
        self.config(cursor='watch')
        self.update()

        value = self.color_value_var.get()
        if value == 'dm':
            idx = 0
        elif value == 'snr':
            idx = 1
        elif value == 'time':
            idx = 2
        else:  # width
            idx = 4

        if self.data.colorProperty != idx:
            self.data.colorProperty = idx
            self.data.draw()

        self.config(cursor='')

    def onColorMap(self):
        """Set the colormap to the specified value and refresh the plots."""
        self.config(cursor='watch')
        self.update()

        name = self.color_map_var.get()

        # Check for inversion
        if self.color_invert_var.get():
            if name == 'gist_gray':
                name = 'gist_yarg'
            else:
                name = '%s_r' % name

        newCM = cm.get_cmap(name)
        if self.data.cmap != newCM:
            self.data.cmap = newCM
            self.data.draw()

        self.config(cursor='')

    def onColorStretch(self):
        """Set the color stretch to the specified scheme and refresh the plots."""
        self.config(cursor='watch')
        self.update()

        stretch = self.color_stretch_var.get()
        if stretch == 'log':
            newN = LogNorm
        elif stretch == 'sqrt':
            newN = SqrtNorm
        elif stretch == 'sqrd':
            newN = SqrdNorm
        elif stretch == 'asinh':
            newN = AsinhNorm
        elif stretch == 'sinh':
            newN = SinhNorm
        elif stretch == 'hist':
            newN = HistEqNorm
        else:  # linear
            newN = Normalize

        if self.data.norm != newN:
            self.data.norm = newN
            self.data.draw()

        self.config(cursor='')

    def onDataSymbol(self):
        """Set the plotting symbol to the specified shape and refresh the plots."""
        self.config(cursor='watch')
        self.update()

        marker = self.plot_symbol_var.get()

        if self.data.plotSymbol != marker:
            self.data.plotSymbol = marker
            self.data.draw()

        self.config(cursor='')

    def onDataSize(self):
        """Set the symbol size coding to the specified quantity and refresh the plots."""
        self.config(cursor='watch')
        self.update()

        value = self.size_map_var.get()
        if value == 'dm':
            idx = 0
        elif value == 'snr':
            idx = 1
        elif value == 'time':
            idx = 2
        else:  # width
            idx = 4

        if self.data.sizeProperty != idx:
            self.data.sizeProperty = idx
            self.data.draw()

        self.config(cursor='')

    def onDataAdjust(self):
        """Bring up the data threshold adjustment dialog window."""
        ThresholdAdjust(self)

    def onDisplayDecimate(self):
        """Toggle scatter plot display decimation on/off."""
        self.config(cursor='watch')
        self.update()

        self.data.fullRes = not self.data.fullRes
        self.data.draw()

        self.config(cursor='')

    def onDisplayDecimateAdjust(self):
        """Bring up the scatter plot display decimation adjustment dialog window."""
        DecimationAdjust(self)

    def onKeyPress(self, event):
        """Handle key presses."""
        keycode = event.keysym
        # Key handling is primarily done in matplotlib event handlers
        pass

    def onSize(self, event):
        """Handle window resize."""
        # Use after_idle to avoid too many redraws
        self.after_idle(self.resizePlots)

    def resizePlots(self):
        """Resize the plots to fit the window."""
        try:
            # Get the current size of the window
            w = self.winfo_width()
            h = self.winfo_height()

            if w < 100 or h < 100:
                return

            # Get toolbar height
            try:
                ht = self.toolbar.winfo_height()
            except:
                ht = 30

            # Come up with new figure sizes in inches
            dpi = self.figure1a.get_dpi()
            newW1 = 1.0 * (w / 3) / dpi
            newH0 = 1.0 * (h / 2) / dpi
            newH1 = 1.0 * (h / 2 - ht) / dpi

            # Set the figure sizes and redraw
            self.figure1a.set_size_inches((newW1, newH0))
            self.figure1a.tight_layout()
            self.canvas1a.draw()

            self.figure1b.set_size_inches((newW1, newH0))
            self.figure1b.tight_layout()
            self.canvas1b.draw()

            self.figure1c.set_size_inches((newW1, newH0))
            self.figure1c.tight_layout()
            self.canvas1c.draw()

            self.figure2.set_size_inches((w / dpi, newH1))
            self.figure2.tight_layout()
            self.canvas2.draw()
        except Exception:
            pass


class ThresholdAdjust(tk.Toplevel):
    """Dialog for adjusting pulse filtering thresholds."""

    def __init__(self, parent):
        super().__init__(parent)
        self.title('Pulse Filtering Adjustment')
        self.transient(parent)
        self.sdf_parent = parent

        self.initUI()
        self.initEvents()

    def initUI(self):
        frame = ttk.Frame(self, padding=10)
        frame.pack(fill=tk.BOTH, expand=True)

        # Min S/N Threshold
        ttk.Label(frame, text='Min. S/N Threshold:').grid(row=0, column=0, sticky=tk.W, padx=5, pady=2)
        self.t_var = tk.StringVar(value='%.1f' % self.sdf_parent.data.dataThreshold[0])
        self.tText = ttk.Entry(frame, textvariable=self.t_var, width=12)
        self.tText.grid(row=0, column=1, padx=5, pady=2)
        ttk.Button(frame, text='-', width=3, command=self.onThresholdDecrease).grid(row=0, column=2, padx=2)
        ttk.Button(frame, text='+', width=3, command=self.onThresholdIncrease).grid(row=0, column=3, padx=2)

        # Max Pulse Width
        ttk.Label(frame, text='Max. Pulse Width [ms]:').grid(row=1, column=0, sticky=tk.W, padx=5, pady=2)
        self.u_var = tk.StringVar(value='%.3f' % self.sdf_parent.data.dataThreshold[2])
        self.uText = ttk.Entry(frame, textvariable=self.u_var, width=12)
        self.uText.grid(row=1, column=1, padx=5, pady=2)
        ttk.Button(frame, text='-', width=3, command=self.onUpperDecrease).grid(row=1, column=2, padx=2)
        ttk.Button(frame, text='+', width=3, command=self.onUpperIncrease).grid(row=1, column=3, padx=2)

        # Min Pulse Width
        ttk.Label(frame, text='Min. Pulse Width [ms]:').grid(row=2, column=0, sticky=tk.W, padx=5, pady=2)
        self.l_var = tk.StringVar(value='%.3f' % self.sdf_parent.data.dataThreshold[1])
        self.lText = ttk.Entry(frame, textvariable=self.l_var, width=12)
        self.lText.grid(row=2, column=1, padx=5, pady=2)
        ttk.Button(frame, text='-', width=3, command=self.onLowerDecrease).grid(row=2, column=2, padx=2)
        ttk.Button(frame, text='+', width=3, command=self.onLowerIncrease).grid(row=2, column=3, padx=2)

        # Separator and OK button
        ttk.Separator(frame, orient=tk.HORIZONTAL).grid(row=3, column=0, columnspan=4, sticky='ew', pady=10)
        ttk.Button(frame, text='Ok', command=self.onOk).grid(row=4, column=3, padx=5, pady=5)

    def initEvents(self):
        self.bind('<Return>', lambda e: self.onKeyPress())

    def onKeyPress(self):
        self.sdf_parent.data.dataThreshold[0] = float(self.t_var.get())
        self.sdf_parent.data.dataThreshold[1] = float(self.l_var.get())
        self.sdf_parent.data.dataThreshold[2] = float(self.u_var.get())
        self.sdf_parent.data.draw(recompute=True)

    def onThresholdDecrease(self):
        self.sdf_parent.data.dataThreshold[0] -= 1
        if self.sdf_parent.data.dataThreshold[0] < 1.0:
            self.sdf_parent.data.dataThreshold[0] = 1.0
        self.t_var.set('%.1f' % self.sdf_parent.data.dataThreshold[0])
        self.sdf_parent.data.draw(recompute=True)

    def onThresholdIncrease(self):
        self.sdf_parent.data.dataThreshold[0] += 1
        self.t_var.set('%.1f' % self.sdf_parent.data.dataThreshold[0])
        self.sdf_parent.data.draw(recompute=True)

    def onUpperDecrease(self):
        self.sdf_parent.data.dataThreshold[2] -= self.sdf_parent.data.dataThreshold[2] * 0.1
        if self.sdf_parent.data.dataThreshold[2] < 0.0:
            self.sdf_parent.data.dataThreshold[2] = 0.0
        self.u_var.set('%.3f' % self.sdf_parent.data.dataThreshold[2])
        self.sdf_parent.data.draw(recompute=True)

    def onUpperIncrease(self):
        self.sdf_parent.data.dataThreshold[2] += self.sdf_parent.data.dataThreshold[2] * 0.1
        self.u_var.set('%.3f' % self.sdf_parent.data.dataThreshold[2])
        self.sdf_parent.data.draw(recompute=True)

    def onLowerDecrease(self):
        self.sdf_parent.data.dataThreshold[1] -= self.sdf_parent.data.dataThreshold[1] * 0.1
        if self.sdf_parent.data.dataThreshold[1] < 0.0:
            self.sdf_parent.data.dataThreshold[1] = 0.0
        self.l_var.set('%.3f' % self.sdf_parent.data.dataThreshold[1])
        self.sdf_parent.data.draw(recompute=True)

    def onLowerIncrease(self):
        self.sdf_parent.data.dataThreshold[1] += self.sdf_parent.data.dataThreshold[1] * 0.1
        self.l_var.set('%.3f' % self.sdf_parent.data.dataThreshold[1])
        self.sdf_parent.data.draw(recompute=True)

    def onOk(self):
        needToRedraw = False

        if self.sdf_parent.data.dataThreshold[0] != float(self.t_var.get()):
            self.sdf_parent.data.dataThreshold[0] = float(self.t_var.get())
            needToRedraw = True
        if self.sdf_parent.data.dataThreshold[1] != float(self.l_var.get()):
            self.sdf_parent.data.dataThreshold[1] = float(self.l_var.get())
            needToRedraw = True
        if self.sdf_parent.data.dataThreshold[2] != float(self.u_var.get()):
            self.sdf_parent.data.dataThreshold[2] = float(self.u_var.get())
            needToRedraw = True

        if needToRedraw:
            self.sdf_parent.data.draw(recompute=True)

        self.destroy()


class DecimationAdjust(tk.Toplevel):
    """Dialog for adjusting display decimation."""

    def __init__(self, parent):
        super().__init__(parent)
        self.title('Display Decimation Adjustment')
        self.transient(parent)
        self.sdf_parent = parent

        self.initUI()
        self.initEvents()

    def initUI(self):
        frame = ttk.Frame(self, padding=10)
        frame.pack(fill=tk.BOTH, expand=True)

        # Max Points to Plot
        ttk.Label(frame, text='Max. Points to Plot:').grid(row=0, column=0, sticky=tk.W, padx=5, pady=2)
        self.u_var = tk.StringVar(value='%i' % self.sdf_parent.data.maxPoints)
        self.uText = ttk.Entry(frame, textvariable=self.u_var, width=12)
        self.uText.grid(row=0, column=1, padx=5, pady=2)
        ttk.Button(frame, text='-', width=3, command=self.onUpperDecrease).grid(row=0, column=2, padx=2)
        ttk.Button(frame, text='+', width=3, command=self.onUpperIncrease).grid(row=0, column=3, padx=2)

        # Separator and OK button
        ttk.Separator(frame, orient=tk.HORIZONTAL).grid(row=1, column=0, columnspan=4, sticky='ew', pady=10)
        ttk.Button(frame, text='Ok', command=self.onOk).grid(row=2, column=3, padx=5, pady=5)

    def initEvents(self):
        self.bind('<Return>', lambda e: self.onKeyPress())

    def onKeyPress(self):
        self.sdf_parent.data.maxPoints = int(self.u_var.get())
        self.sdf_parent.data.draw(recompute=True)

    def onUpperDecrease(self):
        self.sdf_parent.data.maxPoints -= 1000
        if self.sdf_parent.data.maxPoints < 1000:
            self.sdf_parent.data.maxPoints = 1000
        self.u_var.set('%i' % self.sdf_parent.data.maxPoints)

        t0 = time.time()
        self.sdf_parent.data.draw(recompute=True)
        t1 = time.time()
        print("-> Drawing time with %i points is %.3f s" % (self.sdf_parent.data.maxPoints, t1 - t0))

    def onUpperIncrease(self):
        self.sdf_parent.data.maxPoints += 1000
        self.u_var.set('%i' % self.sdf_parent.data.maxPoints)

        t0 = time.time()
        self.sdf_parent.data.draw(recompute=True)
        t1 = time.time()
        print("-> Drawing time with %i points is %.3f s" % (self.sdf_parent.data.maxPoints, t1 - t0))

    def onOk(self):
        needToRedraw = False

        if self.sdf_parent.data.maxPoints != int(self.u_var.get()):
            self.sdf_parent.data.maxPoints = int(self.u_var.get())
            needToRedraw = True

        if needToRedraw:
            self.sdf_parent.data.draw(recompute=True)

        self.destroy()


class SliceDisplay(tk.Toplevel):
    """Window for displaying a time slice of data in a zoomable fashion."""

    def __init__(self, parent, t, dm, width):
        super().__init__(parent)
        self.title('Time Slice')
        self.geometry('400x375')
        self.sdf_parent = parent

        self.t = t
        self.dm = dm
        self.width = width

        self.initUI()
        self.initPlot()

    def initUI(self):
        """Start the user interface."""
        # Status bar
        self.statusbar = ttk.Label(self, text='', relief=tk.SUNKEN, anchor=tk.W)
        self.statusbar.pack(fill=tk.X, side=tk.BOTTOM)

        # Figure frame
        frame = ttk.Frame(self)
        frame.pack(fill=tk.BOTH, expand=True)

        self.figure = Figure()
        self.canvas = FigureCanvasTkAgg(self.figure, master=frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        self.toolbar = NavigationToolbar2Tk(self.canvas, frame)
        self.toolbar.update()

        # Resize binding
        self.bind('<Configure>', self.resizePlots)

    def initPlot(self):
        """Populate the figure/canvas areas with a plot."""
        # Plot the slice
        self.figure.clf()
        self.ax1 = self.figure.gca()

        # Determine what data to plot
        fLow = self.sdf_parent.data.meta.lofreq
        fHigh = self.sdf_parent.data.meta.lofreq + self.sdf_parent.data.meta.BW

        slope = -_D * (1.0 / fLow**2 - 1.0 / fHigh**2)

        sLow, wLow, wHigh = self.sdf_parent.data.dataThreshold
        tCutLow = self.t + (self.sdf_parent.data.dmMax - self.dm) * slope - 1
        tCutHigh = self.t + (self.sdf_parent.data.dmMin - self.dm) * slope + 1

        valid = numpy.where(
            (self.sdf_parent.data.data[:, 2] >= tCutLow) &
            (self.sdf_parent.data.data[:, 2] <= tCutHigh) &
            (self.sdf_parent.data.data[:, 1] >= sLow) &
            (self.sdf_parent.data.data[:, 4] >= wLow) &
            (self.sdf_parent.data.data[:, 4] <= wHigh)
        )[0]

        if len(valid) == 0:
            self.ax1.set_title('Nothing to display')
            return False

        self.subData = self.sdf_parent.data.data[valid, :]

        deltaT = self.subData[:, 2] + (self.subData[:, 0] - self.dm) * slope - self.t
        valid = numpy.where(numpy.abs(deltaT) < 1)[0]
        if len(valid) == 0:
            self.ax1.set_title('Nothing to display')
            return False

        self.subData = self.subData[valid, :]

        m = self.ax1.scatter(
            self.subData[:, 0], self.subData[:, 1],
            c=self.subData[:, self.sdf_parent.data.colorProperty],
            s=self.subData[:, self.sdf_parent.data.sizeProperty] * 5,
            cmap=self.sdf_parent.data.cmap,
            norm=self.sdf_parent.data.norm(*self.sdf_parent.data.limits[self.sdf_parent.data.colorProperty]),
            marker='o', edgecolors='face'
        )

        try:
            cb = self.figure.colorbar(m, use_gridspec=True)
        except:
            if len(self.figure.get_axes()) > 1:
                self.figure.delaxes(self.figure.get_axes()[-1])
            cb = self.figure.colorbar(m)

        if self.sdf_parent.data.colorProperty == 0:
            cb.ax.set_ylabel('DM [pc cm$^{-3}$]')
        elif self.sdf_parent.data.colorProperty == 1:
            cb.ax.set_ylabel('S/N')
        elif self.sdf_parent.data.colorProperty == 2:
            cb.ax.set_ylabel('Elapsed Time [s]')
        else:
            cb.ax.set_ylabel('Width [ms]')

        # Generate the expected S/N vs. DM curve from Cordes & McLaughlin (2003)
        best = numpy.argmax(self.subData[:, 1])
        dm_best = self.subData[best, 0]
        width_best = self.subData[best, 4]

        dms1 = numpy.linspace(self.subData[:, 0].min(), self.subData[:, 0].max(), num=101)
        dms0 = numpy.linspace(self.sdf_parent.data.ax2.get_ylim()[0], dms1.min(), num=101)
        dms2 = numpy.linspace(dms1.max(), self.sdf_parent.data.ax2.get_ylim()[1], num=101)
        dms = numpy.concatenate([dms0, dms1, dms2])
        zeta = 6.91e-3 * (dms - dm_best) * self.sdf_parent.data.meta.BW / width_best / ((fLow + fHigh) / 2.0 / 1000.0)**3
        thr = numpy.sqrt(numpy.pi) / 2 * erf(zeta) / zeta
        thr *= self.subData[:, 1].max()
        self.ax1.plot(dms, thr, linestyle='--', marker='', color='black', label='Cordes & McLaughlin (2003)')
        self.ax1.legend(loc=0)

        self.ax1.axis('auto')
        self.ax1.set_xlim(self.sdf_parent.data.ax2.get_ylim())
        self.ax1.set_ylim((self.subData[:, 1].min() * 0.90, self.subData[:, 1].max() * 1.1))
        self.ax1.set_xlabel('DM [pc cm$^{-3}$]')
        self.ax1.set_ylabel('S/N')
        self.ax1.set_title('Slice through %.1f s, %.3f pc cm$^{-3}$' % (self.t, self.dm))

        self.canvas.draw()
        self.connect()

    def connect(self):
        """Connect to all the events we need."""
        self.cidmotion = self.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)

    def on_motion(self, event):
        """Handle mouse motion events."""
        if event.inaxes:
            clickX = event.xdata
            clickY = event.ydata

            d = (self.subData[:, 0] - clickX)**2 + (self.subData[:, 1] - clickY)**2
            best = numpy.where(
                (d == d.min()) &
                (self.subData[:, 4] >= self.sdf_parent.data.dataThreshold[1]) &
                (self.subData[:, 4] <= self.sdf_parent.data.dataThreshold[2])
            )[0][0]

            dm, snr, width = self.subData[best, 0], self.subData[best, 1], self.subData[best, 4]
            self.statusbar.config(text="DM=%.4f pc cm^-3, S/N=%.1f, width=%.3f ms" % (dm, snr, width))
        else:
            self.statusbar.config(text="")

    def disconnect(self):
        """Disconnect all the stored connection ids."""
        self.figure.canvas.mpl_disconnect(self.cidmotion)

    def resizePlots(self, event=None):
        """Resize plots to fit window."""
        try:
            w = self.winfo_width()
            h = self.winfo_height()
            ht = self.toolbar.winfo_height()

            dpi = self.figure.get_dpi()
            newW = 1.0 * w / dpi
            newH = 1.0 * (h - ht) / dpi
            self.figure.set_size_inches((newW, newH))
            self.figure.tight_layout()
            self.canvas.draw()
        except Exception:
            pass


class WaterfallDisplay(tk.Toplevel):
    """Window for displaying PSRFITS waterfall data."""

    def __init__(self, parent, fitsname, t, dm, width):
        super().__init__(parent)
        self.title('Waterfall')
        self.geometry('600x500')
        self.sdf_parent = parent

        self.fitsname = fitsname
        self.t = t
        self.dm = dm
        self.width = width / 1000.0  # ms -> s

        # Convert time from barycentric to topocentric, if required
        if self.sdf_parent.data.meta.bary:
            if self.sdf_parent.data.bary2topo is not None:
                self.t = self.sdf_parent.data.bary2topo(self.t)

        # Waterfall display settings
        self.cAdjust = None
        self.index = 0
        self.usedB = True
        self.bandpass = True
        self.sweep = True
        self.profile = True
        self.cmap = cm.get_cmap('jet')
        self.norm = Normalize
        self.decFactor = 1

        self.load()

        self.initUI()

        self.render()
        self.draw()

    def initUI(self):
        """Start the user interface."""
        # Status bar
        self.statusbar = ttk.Label(self, text='', relief=tk.SUNKEN, anchor=tk.W)
        self.statusbar.pack(fill=tk.X, side=tk.BOTTOM)

        # Menu bar
        menubar = tk.Menu(self)
        self.config(menu=menubar)

        # Color Menu
        colorMenu = tk.Menu(menubar, tearoff=0)
        colorMenu.add_command(label='Auto-scale Colorbar', command=self.onAutoscale)
        colorMenu.add_command(label='Adjust Contrast', command=self.onColorAdjust)

        ## Color Map submenu
        self.color_map_var = tk.StringVar(value='jet')
        cmapMenu = tk.Menu(colorMenu, tearoff=0)
        cmapMenu.add_radiobutton(label='Paired', variable=self.color_map_var, value='Paired', command=self.onColorMap)
        cmapMenu.add_radiobutton(label='Spectral', variable=self.color_map_var, value='Spectral', command=self.onColorMap)
        cmapMenu.add_radiobutton(label='Bone', variable=self.color_map_var, value='bone', command=self.onColorMap)
        cmapMenu.add_radiobutton(label='Jet', variable=self.color_map_var, value='jet', command=self.onColorMap)
        cmapMenu.add_radiobutton(label='Earth', variable=self.color_map_var, value='gist_earth', command=self.onColorMap)
        cmapMenu.add_radiobutton(label='Heat', variable=self.color_map_var, value='gist_heat', command=self.onColorMap)
        cmapMenu.add_radiobutton(label='NCAR', variable=self.color_map_var, value='gist_ncar', command=self.onColorMap)
        cmapMenu.add_radiobutton(label='Rainbow', variable=self.color_map_var, value='gist_rainbow', command=self.onColorMap)
        cmapMenu.add_radiobutton(label='Stern', variable=self.color_map_var, value='gist_stern', command=self.onColorMap)
        cmapMenu.add_radiobutton(label='Gray', variable=self.color_map_var, value='gist_gray', command=self.onColorMap)
        cmapMenu.add_separator()
        self.color_invert_var = tk.BooleanVar(value=False)
        cmapMenu.add_checkbutton(label='Invert', variable=self.color_invert_var, command=self.onColorMap)
        colorMenu.add_cascade(label='Color Map', menu=cmapMenu)

        ## Color Stretch submenu
        self.color_stretch_var = tk.StringVar(value='linear')
        smapMenu = tk.Menu(colorMenu, tearoff=0)
        smapMenu.add_radiobutton(label='Linear', variable=self.color_stretch_var, value='linear', command=self.onColorStretch)
        smapMenu.add_radiobutton(label='Log', variable=self.color_stretch_var, value='log', command=self.onColorStretch)
        smapMenu.add_radiobutton(label='Square Root', variable=self.color_stretch_var, value='sqrt', command=self.onColorStretch)
        smapMenu.add_radiobutton(label='Squared', variable=self.color_stretch_var, value='sqrd', command=self.onColorStretch)
        smapMenu.add_radiobutton(label='ASinh', variable=self.color_stretch_var, value='asinh', command=self.onColorStretch)
        smapMenu.add_radiobutton(label='Sinh', variable=self.color_stretch_var, value='sinh', command=self.onColorStretch)
        smapMenu.add_radiobutton(label='Histogram Equalization', variable=self.color_stretch_var, value='hist', command=self.onColorStretch)
        colorMenu.add_cascade(label='Color Stretch', menu=smapMenu)

        menubar.add_cascade(label='Color', menu=colorMenu)

        # Data Menu
        dataMenu = tk.Menu(menubar, tearoff=0)
        self.data_product_var = tk.StringVar(value=self.data_products[0])
        for dataProduct in self.data_products:
            dataMenu.add_radiobutton(label=dataProduct, variable=self.data_product_var,
                                     value=dataProduct, command=self.onDataProduct)
        dataMenu.add_separator()
        self.auto_decimation_var = tk.BooleanVar(value=True)
        dataMenu.add_checkbutton(label='Auto Time Decimation', variable=self.auto_decimation_var,
                                 command=self.onAutoDecimation)
        dataMenu.add_command(label='Adjust Time Decimation', command=self.onAdjustDecimation)
        menubar.add_cascade(label='Data', menu=dataMenu)
        self.dataMenu = dataMenu

        # Bandpass Menu
        bandpassMenu = tk.Menu(menubar, tearoff=0)
        self.bandpass_var = tk.BooleanVar(value=True)
        bandpassMenu.add_radiobutton(label='On', variable=self.bandpass_var, value=True, command=self.onBandpass)
        bandpassMenu.add_radiobutton(label='Off', variable=self.bandpass_var, value=False, command=self.onBandpass)
        menubar.add_cascade(label='Bandpass', menu=bandpassMenu)

        # Display Menu
        displayMenu = tk.Menu(menubar, tearoff=0)
        self.sweep_var = tk.BooleanVar(value=True)
        displayMenu.add_checkbutton(label='Show Sweep', variable=self.sweep_var, command=self.onSweep)
        self.profile_var = tk.BooleanVar(value=True)
        displayMenu.add_checkbutton(label='Show Profile', variable=self.profile_var, command=self.onProfile)
        menubar.add_cascade(label='Display', menu=displayMenu)

        # Figure frame
        frame = ttk.Frame(self)
        frame.pack(fill=tk.BOTH, expand=True)

        self.figure = Figure()
        self.canvas = FigureCanvasTkAgg(self.figure, master=frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        self.toolbar = NavigationToolbar2Tk(self.canvas, frame)
        self.toolbar.update()

        # Resize binding
        self.bind('<Configure>', self.resizePlots)

    def load(self):
        """
        Compute and save everything needed for the plot.
        """

        print("Loading PSRFITS metadata...")
        hdulist = astrofits.open(self.fitsname, mode='readonly', memmap=True)

        ## File specifics
        LFFT = hdulist[1].header['NCHAN']
        tInt = hdulist[1].header['TBIN']
        nSubs = hdulist[1].header['NSBLK']
        tSubs = nSubs * tInt
        nPol = hdulist[1].header['NPOL']
        if nPol == 1:
            data_products = ['I',]
            self.usedB = True
        elif nPol == 2:
            if hdulist[0].header['FD_POLN'] == 'CIRC':
                data_products = ['LL', 'RR']
                self.usedB = True
            else:
                data_products = ['XX', 'YY']
                self.usedB = True
        else:
            data_products = ['I', 'Q', 'U', 'V']
            self.usedB = False
        self.data_products = data_products
        nChunks = len(hdulist[1].data)

        ## Frequency information
        freq = hdulist[1].data[0][12] * 1e6
        self.freq = freq

        ## Pulse location
        tSweep = delay(freq, self.dm)
        self.tSweep = tSweep
        tStart = self.t
        tStop = tStart + tSweep.max()
        subIntStart = int(tStart / tSubs) - 1
        subIntStop = int(tStop / tSubs) + 1
        subIntStart = max([0, subIntStart])
        subIntStop = min([subIntStop, nChunks - 1])

        ## Spectra extraction
        print("Extracting event region...")
        samp, tRel, spec, mask = [], [], [], []
        for i in range(subIntStart, subIntStop + 1):
            ### Access the correct subintegration
            subint = hdulist[1].data[i]

            ### Pull out various bits that we need, including:
            ###  * the start time of the subint. - tOff
            ###  * the weight mask, converted to binary - msk
            ###  * the scale and offset values - bscl and bzero
            ###  * the actual data - data
            tOff = subint[1] - subint[0] // 2
            msk = numpy.where(subint[13] >= 0.5, False, True)
            bzero = subint[14]
            bscl = subint[15]
            bzero.shape = (LFFT, nPol)
            bscl.shape = (LFFT, nPol)
            bzero = bzero.T
            bscl = bscl.T
            data = subint[16]
            data.shape = (nSubs, LFFT, nPol)
            data = data.T

            ### Apply the scaling/offset to the data and save the results
            for j in range(nSubs):
                s = i * nSubs + j
                t = subint[1] - self.t + tInt * (j - nSubs // 2)
                d = data[:, :, j] * bscl + bzero
                samp.append(s)
                tRel.append(t)
                spec.append(d)
                mask.append(msk)
        samp = numpy.array(samp)
        self.tRel = numpy.array(tRel)

        self.spec = numpy.ma.array(numpy.array(spec), mask=numpy.array(mask))
        hdulist.close()

        ## Bandpassing
        print("Computing bandpass...")
        try:
            from _helper import FastAxis0Median
            meanSpec = FastAxis0Median(self.spec)
        except ImportError:
            meanSpec = numpy.mean(self.spec, axis=0)

        ### Come up with an appropriate smoothing window (ws) and order (od)
        ws = int(round(self.spec.shape[2] / 10.0))
        ws = min([41, ws])
        if ws % 2 == 0:
            ws += 1
        od = min([9, ws - 2])

        bpm2 = []
        for i in range(self.spec.shape[1]):
            bpm = savitzky_golay(meanSpec[i, :], ws, od, deriv=0)
            bpm = numpy.ma.array(bpm, mask=~numpy.isfinite(bpm))

            if bpm.mean() == 0:
                bpm += 1
            bpm2.append(bpm / bpm.mean())

        ### Apply the bandpass correction
        bpm2 = numpy.array(bpm2)
        self.specBandpass = numpy.ma.array(self.spec.data * 1.0, mask=self.spec.mask)
        try:
            from _helper import FastAxis0Bandpass
            FastAxis0Bandpass(self.specBandpass.data, bpm2.astype(numpy.float32))
        except ImportError:
            for i in range(self.spec.shape[1]):
                self.specBandpass.data[:, i, :] = self.spec.data[:, i, :] / bpm2[i]

        # Downselect to something that centers the pulse
        tPad = 0.2 + self.width
        valid = numpy.where((self.tRel > -tPad) & (self.tRel < (tStop - tStart + tPad)))[0]
        self.tRel = self.tRel[valid]
        self.spec = self.spec[valid, :, :]
        self.specBandpass = self.specBandpass[valid, :, :]

        # Run the incoherent dedispersion on the data
        print("Dedispersing data...")
        self.specD = self.spec * 0
        self.specBandpassD = self.specBandpass * 0
        for i in range(self.specD.shape[1]):
            self.specD[:, i, :] = incoherent(freq, self.spec[:, i, :], tInt, self.dm, boundary='fill', fill_value=numpy.nan)
            self.specBandpassD[:, i, :] = incoherent(freq, self.specBandpass[:, i, :], tInt, self.dm, boundary='fill', fill_value=numpy.nan)

        # Calculate the plot limits
        print("Setting default plot limits...")
        self.limits = [None,] * self.spec.shape[1]
        self.limitsBandpass = [None,] * self.spec.shape[1]

        try:
            from _helper import FastAxis1MinMax
            limits0 = FastAxis1MinMax(self.spec)
            limits1 = FastAxis1MinMax(self.specBandpass, chanMin=self.spec.shape[2] // 10, chanMax=9 * self.spec.shape[2] // 10)
            if self.usedB:
                limits0 = to_dB(limits0)
                limits1 = to_dB(limits1)
            for i in range(self.spec.shape[1]):
                self.limits[i] = list(limits0[i, :])
                self.limitsBandpass[i] = list(limits1[i, :])
        except ImportError:
            toUse = range(self.spec.shape[2] // 10, 9 * self.spec.shape[2] // 10 + 1)
            for i in range(self.spec.shape[1]):
                self.limits[i] = findLimits(self.spec[:, i, :], usedB=self.usedB)
            for i in range(self.spec.shape[1]):
                self.limitsBandpass[i] = findLimits(self.specBandpass[:, i, toUse], usedB=self.usedB)

        # Flip the axis to make the pulsar people happy
        self.spec = self.spec.T
        self.specBandpass = self.specBandpass.T

        # Suggest a time decimation factor
        self.decFactor = max([1, int(round(self.width / 10.0 / self.sdf_parent.data.meta.dt))])

        print("Ready")
        return True

    def render(self):
        """Clear the old figure and connect events."""
        self.figure.clf()
        self.connect()

    def addSubplotAxes(self, fig, ax, rect, axisbg='w'):
        """
        Add a small sub-plot to an existing figure.

        From:
            https://stackoverflow.com/questions/17458580/embedding-small-plots-inside-subplots-in-matplotlib
        """

        box = ax.get_position()
        width = box.width
        height = box.height
        inax_position = ax.transAxes.transform(rect[0:2])
        transFigure = fig.transFigure.inverted()
        infig_position = transFigure.transform(inax_position)
        x = infig_position[0]
        y = infig_position[1]
        width *= rect[2]
        height *= rect[3]
        subax = fig.add_axes([x, y, width, height], facecolor=axisbg)
        x_labelsize = subax.get_xticklabels()[0].get_size()
        y_labelsize = subax.get_yticklabels()[0].get_size()
        x_labelsize *= rect[2] ** 0.5
        y_labelsize *= rect[3] ** 0.5
        subax.xaxis.set_tick_params(labelsize=x_labelsize)
        subax.yaxis.set_tick_params(labelsize=y_labelsize)
        return subax

    def draw(self):
        """
        Populate the figure/canvas areas with a plot.
        """

        # Setup
        tRel = self.tRel.copy()
        tSweep = self.tSweep
        freq = self.freq
        dataProduct = self.data_products[self.index]
        if self.bandpass:
            spec = self.specBandpass[:, self.index, :]
            specD = self.specBandpassD[:, self.index, :]
            limits = self.limitsBandpass[self.index]
        else:
            spec = self.spec[:, self.index, :]
            specD = self.specD[:, self.index, :]
            limits = self.limits[self.index]

        # Decimate
        nKeep = (tRel.size // self.decFactor) * self.decFactor
        tRel = tRel[:nKeep]
        tRel.shape = (nKeep // self.decFactor, self.decFactor)
        tRel = tRel.mean(axis=1)
        spec = spec[:, :nKeep]
        spec.shape = (spec.shape[0], nKeep // self.decFactor, self.decFactor)
        spec = spec.mean(axis=2)

        # Plot Waterfall
        self.figure.clf()
        self.ax1 = self.figure.gca()

        # Get colormap
        cmap_name = self.color_map_var.get()
        if self.color_invert_var.get():
            if cmap_name == 'gist_gray':
                cmap_name = 'gist_yarg'
            else:
                cmap_name = cmap_name + '_r'
        cmap = cm.get_cmap(cmap_name)

        # Get normalization
        stretch = self.color_stretch_var.get()
        if stretch == 'log':
            norm = LogNorm(vmin=limits[0], vmax=limits[1])
        elif stretch == 'sqrt':
            norm = SqrtNorm(vmin=limits[0], vmax=limits[1])
        elif stretch == 'sqrd':
            norm = SqrdNorm(vmin=limits[0], vmax=limits[1])
        elif stretch == 'asinh':
            norm = AsinhNorm(vmin=limits[0], vmax=limits[1])
        elif stretch == 'sinh':
            norm = SinhNorm(vmin=limits[0], vmax=limits[1])
        elif stretch == 'hist':
            norm = HistEqNorm(vmin=limits[0], vmax=limits[1])
        else:
            norm = Normalize(vmin=limits[0], vmax=limits[1])

        if self.usedB:
            m = self.ax1.imshow(to_dB(spec), interpolation='nearest',
                               extent=(tRel[0], tRel[-1], freq[0] / 1e6, freq[-1] / 1e6),
                               origin='lower', cmap=cmap, norm=norm)
            try:
                cb = self.figure.colorbar(m, use_gridspec=True)
            except:
                if len(self.figure.get_axes()) > 1:
                    self.figure.delaxes(self.figure.get_axes()[-1])
                cb = self.figure.colorbar(m)
            cb.ax.set_ylabel('PSD [arb. dB]')
        else:
            m = self.ax1.imshow(spec, interpolation='nearest',
                               extent=(tRel[0], tRel[-1], freq[0] / 1e6, freq[-1] / 1e6),
                               origin='lower', cmap=cmap, norm=norm)
            try:
                cb = self.figure.colorbar(m, use_gridspec=True)
            except:
                if len(self.figure.get_axes()) > 1:
                    self.figure.delaxes(self.figure.get_axes()[-1])
                cb = self.figure.colorbar(m)
            cb.ax.set_ylabel('PSD [arb. lin.]')

        ## Pulse boundary markers
        if self.sweep:
            self.ax1.plot(tSweep - 2 * self.width, freq / 1e6, linestyle='--', color='r')
            self.ax1.plot(tSweep + 2 * self.width, freq / 1e6, linestyle='--', color='r')

        ## Dedispersed profile
        if self.profile:
            prof = specD.sum(axis=0)
            prof = prof[:nKeep]
            prof.shape = (nKeep // self.decFactor, self.decFactor)
            prof = prof.mean(axis=1)
            ax = self.addSubplotAxes(self.figure, self.ax1, [0.7, 0.7, 0.25, 0.25])
            ax.plot(tRel, prof)

        self.ax1.axis('auto')
        self.ax1.set_xlim((tRel[0], tRel[-1]))
        self.ax1.set_ylim((freq[0] / 1e6, freq[-1] / 1e6))
        self.ax1.set_xlabel('Time - %.4f s' % self.t)
        self.ax1.set_ylabel('Frequency [MHz]')
        self.ax1.set_title(dataProduct)

        ## Draw
        self.canvas.draw()

    def connect(self):
        """Connect to events."""
        self.cidmotion = self.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)

    def on_motion(self, event):
        """
        Deal with motion events in the waterfall window. This involves
        setting the status bar with the current x and y coordinates as well
        as the power value.
        """

        if event.inaxes:
            clickX = event.xdata
            clickY = event.ydata

            dataX = numpy.where(numpy.abs(clickY - self.freq / 1e6) == (numpy.abs(clickY - self.freq / 1e6).min()))[0][0]
            dataY = numpy.where(numpy.abs(clickX - self.tRel) == (numpy.abs(clickX - self.tRel).min()))[0][0]

            value = to_dB(self.spec[dataX, 0, dataY])
            self.statusbar.config(text="t=%.6f s, f=%.4f MHz, p=%.2f dB" % (clickX, clickY, value))
        else:
            self.statusbar.config(text="")

    def disconnect(self):
        """Disconnect all the stored connection ids."""
        self.figure.canvas.mpl_disconnect(self.cidmotion)

    def onAutoscale(self):
        """Auto-scale the current data display."""
        self.config(cursor='watch')
        self.update()

        i = self.index
        toUse = numpy.arange(self.spec.shape[0] // 10, 9 * self.spec.shape[0] // 10)

        try:
            from _helper import FastAxis1Percentiles5And99
            if self.bandpass:
                self.limitsBandpass[i] = list(FastAxis1Percentiles5And99(self.specBandpass.T.data, i,
                                              chanMin=self.spec.shape[0] // 10, chanMax=9 * self.spec.shape[0] // 10))
            else:
                self.limits[i] = list(FastAxis1Percentiles5And99(self.spec.T.data, i))
        except ImportError:
            if self.bandpass:
                self.limitsBandpass[i] = numpy.percentile(self.specBandpass[toUse, i, :], (5, 99))
            else:
                self.limits[i] = numpy.percentile(self.spec[:, i, :], (5, 99))

        if self.usedB:
            if self.bandpass:
                self.limitsBandpass[i] = [to_dB(v) for v in self.limitsBandpass[i]]
            else:
                self.limits[i] = [to_dB(v) for v in self.limits[i]]

        self.draw()

        self.config(cursor='')

    def onColorMap(self):
        """Set the colormap to the specified value and refresh the plots."""
        self.config(cursor='watch')
        self.update()

        cmap_name = self.color_map_var.get()
        if self.color_invert_var.get():
            if cmap_name == 'gist_gray':
                cmap_name = 'gist_yarg'
            else:
                cmap_name = cmap_name + '_r'

        newCM = cm.get_cmap(cmap_name)
        if self.cmap != newCM:
            self.cmap = newCM
            self.draw()

        self.config(cursor='')

    def onColorStretch(self):
        """Set the color stretch to the specified scheme and refresh the plots."""
        self.config(cursor='watch')
        self.update()

        stretch = self.color_stretch_var.get()
        if stretch == 'log':
            newN = LogNorm
        elif stretch == 'sqrt':
            newN = SqrtNorm
        elif stretch == 'sqrd':
            newN = SqrdNorm
        elif stretch == 'asinh':
            newN = AsinhNorm
        elif stretch == 'sinh':
            newN = SinhNorm
        elif stretch == 'hist':
            newN = HistEqNorm
        else:
            newN = Normalize

        if self.norm != newN:
            self.norm = newN
            self.draw()

        self.config(cursor='')

    def onColorAdjust(self):
        """Bring up the colorbar adjustment dialog window."""
        WaterfallContrastAdjust(self)

    def onDataProduct(self):
        """Set the data product to display."""
        self.config(cursor='watch')
        self.update()

        newIndex = self.data_products.index(self.data_product_var.get())

        if newIndex != self.index:
            self.index = newIndex
            self.draw()

        self.config(cursor='')

    def onAutoDecimation(self):
        """Automatically adjust the time decimation of the plot."""
        if self.auto_decimation_var.get():
            self.config(cursor='watch')
            self.update()

            # Suggest a time decimation factor
            newDecFactor = max([1, int(round(self.width / 10.0 / self.sdf_parent.data.meta.dt))])

            if newDecFactor != self.decFactor:
                self.decFactor = newDecFactor
                self.draw()

            self.config(cursor='')

    def onAdjustDecimation(self):
        """Manually adjust the time decimation of the plot."""
        WaterfallDecimationAdjust(self)

    def onSweep(self):
        """Toggle the pulse boundaries on/off."""
        self.config(cursor='watch')
        self.update()

        self.sweep = self.sweep_var.get()
        self.draw()

        self.config(cursor='')

    def onProfile(self):
        """Toggle the integrated profile on/off."""
        self.config(cursor='watch')
        self.update()

        self.profile = self.profile_var.get()
        self.draw()

        self.config(cursor='')

    def onBandpass(self):
        """Toggle bandpass correction on/off."""
        self.config(cursor='watch')
        self.update()

        self.bandpass = self.bandpass_var.get()
        self.draw()

        self.config(cursor='')

    def resizePlots(self, event=None):
        """Resize plots to fit window."""
        try:
            w = self.winfo_width()
            h = self.winfo_height()
            ht = self.toolbar.winfo_height()

            dpi = self.figure.get_dpi()
            newW = 1.0 * w / dpi
            newH = 1.0 * (h - ht) / dpi
            self.figure.set_size_inches((newW, newH))
            self.figure.tight_layout()
            self.canvas.draw()
        except Exception:
            pass


class WaterfallContrastAdjust(tk.Toplevel):
    """Dialog for adjusting waterfall contrast."""

    def __init__(self, parent):
        super().__init__(parent)
        self.title('Waterfall Contrast Adjustment')
        self.transient(parent)
        self.sdf_parent = parent

        self.initUI()

        self.sdf_parent.cAdjust = self

    def initUI(self):
        frame = ttk.Frame(self, padding=10)
        frame.pack(fill=tk.BOTH, expand=True)

        pol = self.sdf_parent.data_products[self.sdf_parent.index]
        if self.sdf_parent.bandpass:
            ttk.Label(frame, text='%s - Bandpass' % pol).grid(row=0, column=0, columnspan=4, pady=5)
        else:
            ttk.Label(frame, text='%s' % pol).grid(row=0, column=0, columnspan=4, pady=5)

        # Upper Limit
        ttk.Label(frame, text='Upper Limit:').grid(row=1, column=0, sticky=tk.W, padx=5, pady=2)
        if self.sdf_parent.bandpass:
            self.u_var = tk.StringVar(value='%.1f' % self.sdf_parent.limitsBandpass[self.sdf_parent.index][1])
        else:
            self.u_var = tk.StringVar(value='%.1f' % self.sdf_parent.limits[self.sdf_parent.index][1])
        self.uText = ttk.Entry(frame, textvariable=self.u_var, width=12)
        self.uText.grid(row=1, column=1, padx=5, pady=2)
        ttk.Button(frame, text='-', width=3, command=self.onUpperDecrease).grid(row=1, column=2, padx=2)
        ttk.Button(frame, text='+', width=3, command=self.onUpperIncrease).grid(row=1, column=3, padx=2)

        # Lower Limit
        ttk.Label(frame, text='Lower Limit:').grid(row=2, column=0, sticky=tk.W, padx=5, pady=2)
        if self.sdf_parent.bandpass:
            self.l_var = tk.StringVar(value='%.1f' % self.sdf_parent.limitsBandpass[self.sdf_parent.index][0])
        else:
            self.l_var = tk.StringVar(value='%.1f' % self.sdf_parent.limits[self.sdf_parent.index][0])
        self.lText = ttk.Entry(frame, textvariable=self.l_var, width=12)
        self.lText.grid(row=2, column=1, padx=5, pady=2)
        ttk.Button(frame, text='-', width=3, command=self.onLowerDecrease).grid(row=2, column=2, padx=2)
        ttk.Button(frame, text='+', width=3, command=self.onLowerIncrease).grid(row=2, column=3, padx=2)

        # Range (read-only)
        ttk.Label(frame, text='Range:').grid(row=3, column=0, sticky=tk.W, padx=5, pady=2)
        self.r_var = tk.StringVar(value='%.1f' % self.__getRange())
        self.rText = ttk.Entry(frame, textvariable=self.r_var, width=12, state='readonly')
        self.rText.grid(row=3, column=1, padx=5, pady=2)

        # Separator and OK button
        ttk.Separator(frame, orient=tk.HORIZONTAL).grid(row=4, column=0, columnspan=4, sticky='ew', pady=10)
        ttk.Button(frame, text='Ok', command=self.onOk).grid(row=5, column=3, padx=5, pady=5)

        self.bind('<Return>', lambda e: self.onKeyPress())

    def __getRange(self):
        index = self.sdf_parent.index
        if self.sdf_parent.bandpass:
            return self.sdf_parent.limitsBandpass[index][1] - self.sdf_parent.limitsBandpass[index][0]
        else:
            return self.sdf_parent.limits[index][1] - self.sdf_parent.limits[index][0]

    def __getIncrement(self):
        return 0.1 * self.__getRange()

    def onKeyPress(self):
        index = self.sdf_parent.index
        if self.sdf_parent.bandpass:
            self.sdf_parent.limitsBandpass[index][0] = float(self.l_var.get())
            self.sdf_parent.limitsBandpass[index][1] = float(self.u_var.get())
        else:
            self.sdf_parent.limits[index][0] = float(self.l_var.get())
            self.sdf_parent.limits[index][1] = float(self.u_var.get())
        self.r_var.set('%.1f' % self.__getRange())
        self.sdf_parent.draw()

    def onUpperDecrease(self):
        index = self.sdf_parent.index
        if self.sdf_parent.bandpass:
            self.sdf_parent.limitsBandpass[index][1] -= self.__getIncrement()
            self.u_var.set('%.1f' % self.sdf_parent.limitsBandpass[index][1])
        else:
            self.sdf_parent.limits[index][1] -= self.__getIncrement()
            self.u_var.set('%.1f' % self.sdf_parent.limits[index][1])
        self.r_var.set('%.1f' % self.__getRange())
        self.sdf_parent.draw()

    def onUpperIncrease(self):
        index = self.sdf_parent.index
        if self.sdf_parent.bandpass:
            self.sdf_parent.limitsBandpass[index][1] += self.__getIncrement()
            self.u_var.set('%.1f' % self.sdf_parent.limitsBandpass[index][1])
        else:
            self.sdf_parent.limits[index][1] += self.__getIncrement()
            self.u_var.set('%.1f' % self.sdf_parent.limits[index][1])
        self.r_var.set('%.1f' % self.__getRange())
        self.sdf_parent.draw()

    def onLowerDecrease(self):
        index = self.sdf_parent.index
        if self.sdf_parent.bandpass:
            self.sdf_parent.limitsBandpass[index][0] -= self.__getIncrement()
            self.l_var.set('%.1f' % self.sdf_parent.limitsBandpass[index][0])
        else:
            self.sdf_parent.limits[index][0] -= self.__getIncrement()
            self.l_var.set('%.1f' % self.sdf_parent.limits[index][0])
        self.r_var.set('%.1f' % self.__getRange())
        self.sdf_parent.draw()

    def onLowerIncrease(self):
        index = self.sdf_parent.index
        if self.sdf_parent.bandpass:
            self.sdf_parent.limitsBandpass[index][0] += self.__getIncrement()
            self.l_var.set('%.1f' % self.sdf_parent.limitsBandpass[index][0])
        else:
            self.sdf_parent.limits[index][0] += self.__getIncrement()
            self.l_var.set('%.1f' % self.sdf_parent.limits[index][0])
        self.r_var.set('%.1f' % self.__getRange())
        self.sdf_parent.draw()

    def onOk(self):
        needToRedraw = False
        index = self.sdf_parent.index

        if self.sdf_parent.bandpass:
            if float(self.l_var.get()) != self.sdf_parent.limitsBandpass[index][0]:
                self.sdf_parent.limitsBandpass[index][0] = float(self.l_var.get())
                needToRedraw = True
            if float(self.u_var.get()) != self.sdf_parent.limitsBandpass[index][1]:
                self.sdf_parent.limitsBandpass[index][1] = float(self.u_var.get())
                needToRedraw = True
        else:
            if float(self.l_var.get()) != self.sdf_parent.limits[index][0]:
                self.sdf_parent.limits[index][0] = float(self.l_var.get())
                needToRedraw = True
            if float(self.u_var.get()) != self.sdf_parent.limits[index][1]:
                self.sdf_parent.limits[index][1] = float(self.u_var.get())
                needToRedraw = True

        if needToRedraw:
            self.sdf_parent.draw()

        self.sdf_parent.cAdjust = None
        self.destroy()


class WaterfallDecimationAdjust(tk.Toplevel):
    """Dialog for adjusting waterfall time decimation."""

    def __init__(self, parent):
        super().__init__(parent)
        self.title('Waterfall Decimation Adjustment')
        self.transient(parent)
        self.sdf_parent = parent

        self.initUI()

    def initUI(self):
        frame = ttk.Frame(self, padding=10)
        frame.pack(fill=tk.BOTH, expand=True)

        # Time Decimation
        ttk.Label(frame, text='Time Decimation:').grid(row=0, column=0, sticky=tk.W, padx=5, pady=2)
        self.u_var = tk.StringVar(value='%i' % self.sdf_parent.decFactor)
        self.uText = ttk.Entry(frame, textvariable=self.u_var, width=12)
        self.uText.grid(row=0, column=1, padx=5, pady=2)
        ttk.Button(frame, text='-', width=3, command=self.onDecrease).grid(row=0, column=2, padx=2)
        ttk.Button(frame, text='+', width=3, command=self.onIncrease).grid(row=0, column=3, padx=2)

        # Bins/Pulse (read-only)
        ttk.Label(frame, text='Bins/Pulse:').grid(row=1, column=0, sticky=tk.W, padx=5, pady=2)
        self.r_var = tk.StringVar(value='%.2f' % self.__getBinsPerPulse())
        self.rText = ttk.Entry(frame, textvariable=self.r_var, width=12, state='readonly')
        self.rText.grid(row=1, column=1, padx=5, pady=2)

        # Separator and OK button
        ttk.Separator(frame, orient=tk.HORIZONTAL).grid(row=2, column=0, columnspan=4, sticky='ew', pady=10)
        ttk.Button(frame, text='Ok', command=self.onOk).grid(row=3, column=3, padx=5, pady=5)

        self.bind('<Return>', lambda e: self.onKeyPress())

    def __getBinsPerPulse(self):
        return self.sdf_parent.width / (self.sdf_parent.decFactor * self.sdf_parent.sdf_parent.data.meta.dt)

    def onKeyPress(self):
        self.sdf_parent.decFactor = int(self.u_var.get())
        self.r_var.set('%.2f' % self.__getBinsPerPulse())
        self.sdf_parent.draw()

    def onDecrease(self):
        self.sdf_parent.decFactor -= 1
        if self.sdf_parent.decFactor < 1:
            self.sdf_parent.decFactor = 1
        self.u_var.set('%i' % self.sdf_parent.decFactor)
        self.r_var.set('%.2f' % self.__getBinsPerPulse())
        self.sdf_parent.draw()

    def onIncrease(self):
        self.sdf_parent.decFactor += 1
        self.u_var.set('%i' % self.sdf_parent.decFactor)
        self.r_var.set('%.2f' % self.__getBinsPerPulse())
        self.sdf_parent.draw()

    def onOk(self):
        needToRedraw = False

        if self.sdf_parent.decFactor != int(self.u_var.get()):
            self.sdf_parent.decFactor = int(self.u_var.get())
            needToRedraw = True

        if needToRedraw:
            self.sdf_parent.draw()

        self.destroy()


class HelpWindow(tk.Toplevel):
    """Help window displaying formatted help content."""

    def __init__(self, parent):
        super().__init__(parent)
        self.title('plotSinglePulse Handbook')
        self.geometry('570x400')

        self.anchors = {}
        self.initUI()

    def initUI(self):
        # Create scrollable text widget
        frame = ttk.Frame(self)
        frame.pack(fill=tk.BOTH, expand=True)

        self.text = tk.Text(frame, wrap=tk.WORD, padx=10, pady=10)
        scrollbar = ttk.Scrollbar(frame, orient=tk.VERTICAL, command=self.text.yview)
        self.text.configure(yscrollcommand=scrollbar.set)

        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        self.text.pack(fill=tk.BOTH, expand=True)

        self._configure_tags()
        self._load_help_content()
        self.text.config(state='disabled')

        # Status bar
        self.statusbar = ttk.Label(self, text='', relief=tk.SUNKEN, anchor=tk.W)
        self.statusbar.pack(fill=tk.X, side=tk.BOTTOM)

    def _configure_tags(self):
        default_font = tkfont.nametofont('TkDefaultFont')
        base_size = default_font.cget('size')
        if base_size < 0:
            base_size = 12
        base_family = default_font.cget('family')

        self.text.tag_configure('h4', font=(base_family, base_size + 4, 'bold'), spacing3=5)
        self.text.tag_configure('h6', font=(base_family, base_size + 2, 'bold'), spacing1=10, spacing3=3)
        self.text.tag_configure('bold', font=(base_family, base_size, 'bold'))
        self.text.tag_configure('link', foreground='blue', underline=True)
        self.text.tag_bind('link', '<Enter>', lambda e: self.text.config(cursor='hand2'))
        self.text.tag_bind('link', '<Leave>', lambda e: self.text.config(cursor=''))
        self.text.tag_configure('listitem', lmargin1=20, lmargin2=35)

    def _load_help_content(self):
        """Load help content as formatted text."""
        t = self.text

        # Table of Contents
        self.anchors['top'] = t.index('end')
        t.insert('end', 'Table of Contents\n', 'h4')
        self._insert_link('Introduction', 'intro')
        self._insert_link('Window Layout', 'layout')
        self._insert_link('Usage', 'usage')
        self._insert_link('Mouse Interaction', 'mouse')
        self._insert_link('Keyboard Interaction', 'keyboard')
        t.insert('end', '\n')

        # Introduction
        self.anchors['intro'] = t.index('end')
        t.insert('end', 'Introduction\n', 'h6')
        t.insert('end', 'plotSinglePulse is a graphical interface for working with single pulse search '
                       'results from PRESTO. This script allows the user to plot the candidates in an '
                       'interactive fashion similar to the postscript plots generated by single_pulse_search.py '
                       'as well as to plot the associated waterfall data if a PSRFITS file is provided.\n')
        self._insert_back_link()

        # Window Layout
        self.anchors['layout'] = t.index('end')
        t.insert('end', 'Window Layout\n', 'h6')
        t.insert('end', 'The plotSinglePulse window is broken into two vertical sections. The top section shows, '
                       'from left to right, the histogram of signal-to-noise ratio (S/N) for all pulses displayed, '
                       'the histogram of dispersion measure (DM) for all pulses displayed, and a scatter plot of '
                       'DM vs. S/N. The lower section shows a scatter plot of time vs. DM with user-configurable '
                       'colors and symbols.\n\n')
        t.insert('end', 'Note: ', 'bold')
        t.insert('end', 'Since there may be many candidates the upper left and lower plot windows use internal '
                       'decimation to display the data. This default behavior limits the number of points plotted '
                       'in each window and the decimation can be adjusted via the Data menu.\n')
        self._insert_back_link()

        # Usage
        self.anchors['usage'] = t.index('end')
        t.insert('end', 'Usage\n', 'h6')
        t.insert('end', 'After the single pulse candidates are loaded there are a variety of menu, mouse, '
                       'and keyboard commands that can be used to interact with the data. The menus are:\n')
        t.insert('end', '\u2022 Color - Adjust the parameter used for color coding, color stretch, map, and transfer function\n', 'listitem')
        t.insert('end', '\u2022 Data - Change the symbol coding, size, and pulse selection thresholds\n', 'listitem')
        t.insert('end', '\u2022 Display - Toggle display decimation on/off and change the decimation factors\n', 'listitem')
        t.insert('end', '\u2022 Help - Show this help message\n', 'listitem')
        self._insert_back_link()

        # Mouse Interaction
        self.anchors['mouse'] = t.index('end')
        t.insert('end', 'Mouse Interaction\n', 'h6')
        t.insert('end', 'The mouse can be used in any of the plotting panels. In the upper left and center sections:\n')
        t.insert('end', '\u2022 Left Click - Nothing\n', 'listitem')
        t.insert('end', '\u2022 Middle Click - Select pulses in the current bin for highlighting in the lower panel\n', 'listitem')
        t.insert('end', '\u2022 Right Click - Clear the last selection made with the middle click\n', 'listitem')
        t.insert('end', '\nIn the lower section:\n')
        t.insert('end', '\u2022 Left Click - Move the cross hairs around the screen\n', 'listitem')
        t.insert('end', '\u2022 Middle Click - Unmask the closest pulse candidate to the cross hairs\n', 'listitem')
        t.insert('end', '\u2022 Right Click - Mask the closest pulse candidate to the cross hairs\n', 'listitem')
        t.insert('end', '\nThe lower section can also be zoomed by selecting the zoom tool from the toolbar. ')
        t.insert('end', 'Note: ', 'bold')
        t.insert('end', 'While the zoom tool is active the standard mouse functions in this window are disabled.\n')
        self._insert_back_link()

        # Keyboard Interaction
        self.anchors['keyboard'] = t.index('end')
        t.insert('end', 'Keyboard Interaction\n', 'h6')
        t.insert('end', 'The keyboard can also be used to interact with the lower plotting panel. The keyboard commands are:\n')
        t.insert('end', '\u2022 p - print the information about the underlying pulse\n', 'listitem')
        t.insert('end', '\u2022 e - export the current pulses to a file\n', 'listitem')
        t.insert('end', '\u2022 s - display a DM time slice\n', 'listitem')
        t.insert('end', '\u2022 w - display the PSRFITS waterfall for a pulse\n', 'listitem')
        t.insert('end', '\u2022 u - unmask pulses in a region of time\n', 'listitem')
        t.insert('end', '\u2022 m - mask pulses in a region of time\n', 'listitem')
        t.insert('end', '\u2022 y - unmask pulses in a region of time/DM\n', 'listitem')
        t.insert('end', '\u2022 n - mask pulses in a region of time/DM\n', 'listitem')
        t.insert('end', '\u2022 t - unmask pulses in a region of DM\n', 'listitem')
        t.insert('end', '\u2022 b - mask pulses in a region of DM\n', 'listitem')
        t.insert('end', '\u2022 h - print help message\n', 'listitem')
        t.insert('end', '\n')
        t.insert('end', 'Note: ', 'bold')
        t.insert('end', 'In order to interact with the lower section you will need to click in the axes with the left mouse button.\n')
        self._insert_back_link()

    def _insert_link(self, text, anchor):
        """Insert a clickable link."""
        tag_name = f'link_{anchor}'
        self.text.tag_configure(tag_name, foreground='blue', underline=True)
        self.text.tag_bind(tag_name, '<Button-1>', lambda e, a=anchor: self._scroll_to_anchor(a))
        self.text.tag_bind(tag_name, '<Enter>', lambda e: self.text.config(cursor='hand2'))
        self.text.tag_bind(tag_name, '<Leave>', lambda e: self.text.config(cursor=''))
        self.text.insert('end', f'\u2022 {text}\n', ('listitem', tag_name))

    def _insert_back_link(self):
        """Insert a 'Back to Top' link."""
        tag_name = 'link_back_to_top'
        if not self.text.tag_names().__contains__(tag_name):
            self.text.tag_configure(tag_name, foreground='blue', underline=True)
            self.text.tag_bind(tag_name, '<Button-1>', lambda e: self._scroll_to_anchor('top'))
            self.text.tag_bind(tag_name, '<Enter>', lambda e: self.text.config(cursor='hand2'))
            self.text.tag_bind(tag_name, '<Leave>', lambda e: self.text.config(cursor=''))
        self.text.insert('end', 'Top\n\n', tag_name)

    def _scroll_to_anchor(self, anchor):
        """Scroll to a named anchor."""
        if anchor in self.anchors:
            self.text.see(self.anchors[anchor])


def main(args):
    # Turn off all NumPy warnings to keep stdout clean
    numpy.seterr(all='ignore', invalid='ignore', divide='ignore')

    # Parse the command line options
    args.time_range = [float(v) for v in args.time_range.split(',')]
    args.dm_range = [float(v) for v in args.dm_range.split(',')]
    args.width_range = [float(v) for v in args.width_range.split(',')]

    # Check for the _helper module
    try:
        import _helper
    except ImportError:
        print("WARNING: _helper.so not found, consider building it with 'make'")

    # Go!
    app = MainWindow()
    app.data = SinglePulse_GUI(app)
    app.render()

    if args.filename is not None:
        # If there is a filename on the command line, load it
        app.filenames = args.filename
        app.data.loadData(args.filename, threshold=args.threshold,
                         timeRange=args.time_range, dmRange=args.dm_range,
                         widthRange=args.width_range, fitsname=args.fitsname)
        app.data.render()
        app.data.draw()

    app.mainloop()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='read in a collection of .singlepulse files and plot them interactively',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('filename', type=str, nargs='*', default=None,
                        help='filename to display')
    parser.add_argument('-t', '--threshold', type=aph.positive_float, default=5.0,
                        help='minimum pulsar threshold to display in sigma')
    parser.add_argument('-r', '--time-range', type=str, default='0,inf',
                        help='comma separated list of the relative time range in seconds to load')
    parser.add_argument('-d', '--dm-range', type=str, default='0,inf',
                        help='comma separated list of the DM range in pc cm^-3 to load')
    parser.add_argument('-w', '--width-range', type=str, default='0,inf',
                        help='comma separated list of the pulse width range in ms to load')
    parser.add_argument('-f', '--fitsname', type=str,
                        help='optional PSRFITS file to use for waterfall plots')
    args = parser.parse_args()
    main(args)
