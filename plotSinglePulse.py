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

import wx
import wx.html
import matplotlib
matplotlib.use('WXAgg')
matplotlib.interactive(True)

from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg, FigureCanvasWxAgg
from matplotlib.colors import Normalize
from matplotlib.collections import CircleCollection
from matplotlib import cm
from matplotlib.figure import Figure

__version__ = "0.1"
__author__ = "Jayce Dowell"


# Deal with the different wxPython versions
if 'phoenix' in wx.PlatformInfo:
    AppendMenuItem = lambda x, y: x.Append(y)
    AppendMenuMenu = lambda *args, **kwds: args[0].Append(*args[1:], **kwds)
else:
    AppendMenuItem = lambda x, y: x.AppendItem(y)
    AppendMenuMenu = lambda *args, **kwds: args[0].AppendMenu(*args[1:], **kwds)


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


class RefreshAwareToolbar(NavigationToolbar2WxAgg):
    """
    Sub-class of NavigationToolbar2WxAgg that includes a reference to a 
    callback function that is called when the plot window is changed by the 
    various buttons.  This sub-class also excludes the plot plot button 
    courtesy of:
    
    https://stackoverflow.com/questions/12695678/how-to-modify-the-navigation-toolbar-easily-in-a-matplotlib-figure-window
    """
    
    # only display the buttons we need
    toolitems = [t for t in NavigationToolbar2WxAgg.toolitems if t[0] in ('Home', 'Forward', 'Back', 'Pan', 'Zoom', 'Save')]
    
    def __init__(self, plotCanvas, refreshCallback=None):
        NavigationToolbar2WxAgg.__init__(self, plotCanvas)
        self.refreshCallback = refreshCallback
        POSITION_OF_CONFIGURE_SUBPLOTS_BTN = 6
        try:
            self.DeleteToolByPos(POSITION_OF_CONFIGURE_SUBPLOTS_BTN)
        except wx.wxAssertionError:
            pass
        
    def home(self, *args):
        NavigationToolbar2WxAgg.home(self, *args)
        
        if self.refreshCallback is not None:
            self.refreshCallback()
            
    def forward(self, *args):
        NavigationToolbar2WxAgg.forward(self, *args)
        
        if self.refreshCallback is not None:
            self.refreshCallback()
            
    def back(self, *args):
        NavigationToolbar2WxAgg.back(self, *args)
        
        if self.refreshCallback is not None:
            self.refreshCallback()
            
    def release_zoom(self, *args):
        NavigationToolbar2WxAgg.release_zoom(self, *args)
        
        if self.refreshCallback is not None:
            self.refreshCallback()
            
    def release_pan(self, *args):
        NavigationToolbar2WxAgg.release_pan(self, *args)
        
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
        
            self.frame.statusbar.SetStatusText("%.3f s, %.3f pc cm^-3" % (clickX, clickY))
        else:
            self.frame.statusbar.SetStatusText("")
            
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


ID_OPEN    = 10
ID_QUIT    = 11

ID_COLOR_VALUE_DM = 2100
ID_COLOR_VALUE_SNR = 2101
ID_COLOR_VALUE_TIME = 2102
ID_COLOR_VALUE_WIDTH = 2103
ID_COLOR_MAP_PAIRED = 2200
ID_COLOR_MAP_SPECTRAL = 2201
ID_COLOR_MAP_BONE = 2202
ID_COLOR_MAP_JET = 2203
ID_COLOR_MAP_EARTH = 2204
ID_COLOR_MAP_HEAT = 2205
ID_COLOR_MAP_NCAR = 2206
ID_COLOR_MAP_RAINBOW = 2207
ID_COLOR_MAP_STERN = 2208
ID_COLOR_MAP_GRAY = 2209
ID_COLOR_INVERT = 23
ID_COLOR_STRETCH_LINEAR = 2400
ID_COLOR_STRETCH_LOG = 2401
ID_COLOR_STRETCH_SQRT = 2402
ID_COLOR_STRETCH_SQRD = 2403
ID_COLOR_STRETCH_ASINH = 2404
ID_COLOR_STRETCH_SINH = 2405
ID_COLOR_STRETCH_HIST = 2406

ID_DATA_SYMBOL_CIRCLE = 3000
ID_DATA_SYMBOL_SQUARE = 3001
ID_DATA_SYMBOL_DIAMOND = 3002
ID_DATA_SYMBOL_HEXAGON = 3003
ID_DATA_SYMBOL_PLUS = 3004
ID_DATA_SIZE_DM = 3100
ID_DATA_SIZE_SNR = 3101
ID_DATA_SIZE_TIME = 3102
ID_DATA_SIZE_WIDTH = 3103
ID_DATA_ADJUST = 32

ID_DISPLAY_DECIMATE = 40
ID_DISPLAY_DECIMATE_ADJUST = 41

ID_HELP = 70
ID_ABOUT = 71

class MainWindow(wx.Frame):
    def __init__(self, parent, id):
        self.dirname = ''
        self.filenames = []
        self.data = None
        self.examineFileButton = None
        self.examineWindow = None
        
        wx.Frame.__init__(self, parent, id, title="Single Pulse Viewer", size=(1000,600))
        
    def render(self):
        self.initUI()
        self.initEvents()
        self.Show()
        self.SetClientSize((1000,600))
        
    def initUI(self):
        self.statusbar = self.CreateStatusBar() # A Statusbar in the bottom of the window
        
        font = wx.SystemSettings.GetFont(wx.SYS_SYSTEM_FONT)
        font.SetPointSize(10)
        
        menuBar = wx.MenuBar()
        
        fileMenu = wx.Menu()
        colorMenu = wx.Menu()
        dataMenu = wx.Menu()
        displayMenu = wx.Menu()
        helpMenu = wx.Menu()
        
        ## File Menu
        open = wx.MenuItem(fileMenu, ID_OPEN, "&Open")
        AppendMenuItem(fileMenu, open)
        fileMenu.AppendSeparator()
        exit = wx.MenuItem(fileMenu, ID_QUIT, "E&xit")
        AppendMenuItem(fileMenu, exit)
        
        ## Color Menu
        vmap = wx.Menu()
        vmap.AppendRadioItem(ID_COLOR_VALUE_DM, '&DM')
        vmap.AppendRadioItem(ID_COLOR_VALUE_SNR, '&S/N')
        vmap.AppendRadioItem(ID_COLOR_VALUE_TIME, '&Time')
        vmap.AppendRadioItem(ID_COLOR_VALUE_WIDTH, '&Width')
        AppendMenuMenu(colorMenu, -1, 'Color &Mapping', vmap)
        vmap.Check(ID_COLOR_VALUE_WIDTH, True)
        self.vmapMenu = vmap
        cmap = wx.Menu()
        cmap.AppendRadioItem(ID_COLOR_MAP_PAIRED, '&Paired')
        cmap.AppendRadioItem(ID_COLOR_MAP_SPECTRAL, "&Spectral")
        cmap.AppendRadioItem(ID_COLOR_MAP_BONE, '&Bone')
        cmap.AppendRadioItem(ID_COLOR_MAP_JET, '&Jet')
        cmap.AppendRadioItem(ID_COLOR_MAP_EARTH, '&Earth')
        cmap.AppendRadioItem(ID_COLOR_MAP_HEAT, '&Heat')
        cmap.AppendRadioItem(ID_COLOR_MAP_NCAR, '&NCAR')
        cmap.AppendRadioItem(ID_COLOR_MAP_RAINBOW, '&Rainbow')
        cmap.AppendRadioItem(ID_COLOR_MAP_STERN, 'S&tern')
        cmap.AppendRadioItem(ID_COLOR_MAP_GRAY, '&Gray')
        cmap.AppendSeparator()
        cmap.AppendCheckItem(ID_COLOR_INVERT, 'In&vert')
        AppendMenuMenu(colorMenu, -1, 'Color &Map', cmap)
        cmap.Check(ID_COLOR_MAP_JET, True)
        self.cmapMenu = cmap
        smap = wx.Menu()
        smap.AppendRadioItem(ID_COLOR_STRETCH_LINEAR, '&Linear')
        smap.AppendRadioItem(ID_COLOR_STRETCH_LOG, 'Lo&g')
        smap.AppendRadioItem(ID_COLOR_STRETCH_SQRT, 'Square &Root')
        smap.AppendRadioItem(ID_COLOR_STRETCH_SQRD, '&Squared')
        smap.AppendRadioItem(ID_COLOR_STRETCH_ASINH, '&ASinh')
        smap.AppendRadioItem(ID_COLOR_STRETCH_SINH, '&Sinh')
        smap.AppendRadioItem(ID_COLOR_STRETCH_HIST, '&Histogram Equalization')
        AppendMenuMenu(colorMenu, -1, 'Color &Stretch', smap)
        smap.Check(ID_COLOR_STRETCH_LINEAR, True)
        self.smapMenu = smap
        
        ## Data Menu
        mmap = wx.Menu()
        mmap.AppendRadioItem(ID_DATA_SYMBOL_CIRCLE, '&Circle')
        mmap.AppendRadioItem(ID_DATA_SYMBOL_SQUARE, '&Square')
        mmap.AppendRadioItem(ID_DATA_SYMBOL_DIAMOND, '&Diamond')
        mmap.AppendRadioItem(ID_DATA_SYMBOL_HEXAGON, '&Hexagon')
        mmap.AppendRadioItem(ID_DATA_SYMBOL_PLUS, '&Plus Sign')
        dataMenu.AppendMenu(-1, '&Plot Symbol', mmap)
        mmap.Check(ID_DATA_SYMBOL_CIRCLE, True)
        self.mmapMenu = mmap
        amap = wx.Menu()
        amap.AppendRadioItem(ID_DATA_SIZE_DM, '&DM')
        amap.AppendRadioItem(ID_DATA_SIZE_SNR, '&S/N')
        amap.AppendRadioItem(ID_DATA_SIZE_TIME, '&Time')
        amap.AppendRadioItem(ID_DATA_SIZE_WIDTH, '&Width')
        dataMenu.AppendMenu(-1, 'Size &Mapping', amap)
        amap.Check(ID_DATA_SIZE_SNR, True)
        self.amapMenu = amap
        tadj = wx.MenuItem(dataMenu, ID_DATA_ADJUST, 'Adjust &Thresholds')
        AppendMenuItem(dataMenu, tadj)
        
        ## Display Menu
        displayMenu.AppendCheckItem(ID_DISPLAY_DECIMATE, '&Decimation')
        dadj = wx.MenuItem(dataMenu, ID_DISPLAY_DECIMATE_ADJUST, '&Decimation Adjust')
        AppendMenuItem(displayMenu, dadj)
        self.dadj = dadj
        if not self.data.fullRes:
            displayMenu.Check(ID_DISPLAY_DECIMATE, True)
            dadj.Enable(True)
        else:
            dadj.Enable(False)
            
        ## Help menu items
        help = wx.MenuItem(helpMenu, ID_HELP, 'plotSinglePulse Handbook\tF1')
        AppendMenuItem(helpMenu, help)
        helpMenu.AppendSeparator()
        about = wx.MenuItem(helpMenu, ID_ABOUT, '&About')
        AppendMenuItem(helpMenu, about)
        
        # Creating the menubar.
        menuBar.Append(fileMenu,"&File") # Adding the "filemenu" to the MenuBar
        menuBar.Append(colorMenu, "&Color")
        menuBar.Append(dataMenu, "&Data")
        menuBar.Append(displayMenu, "D&isplay")
        menuBar.Append(helpMenu, "&Help")
        self.SetMenuBar(menuBar)  # Adding the MenuBar to the Frame content.
        
        # Save the various menus so that we can disable them if need be
        self.fileMenu = fileMenu
        self.colorMenu = colorMenu
        self.dataMenu = dataMenu
        self.displayMenu = displayMenu
        self.helpMenu = helpMenu
        
        vbox = wx.BoxSizer(wx.VERTICAL)
        
        # Add SNR histogram plot
        panel1 = wx.Panel(self, -1)
        hbox1 = wx.BoxSizer(wx.HORIZONTAL)
        self.figure1a = Figure(figsize=(2,2))
        self.canvas1a = FigureCanvasWxAgg(panel1, -1, self.figure1a)
        hbox1.Add(self.canvas1a, 1, wx.EXPAND)
        
        # Add a DM histogram plot
        self.figure1b = Figure(figsize=(2,2))
        self.canvas1b = FigureCanvasWxAgg(panel1, -1, self.figure1b)
        hbox1.Add(self.canvas1b, 1, wx.EXPAND)
        
        # Add DM vs SNR plot
        self.figure1c = Figure(figsize=(2,2))
        self.canvas1c = FigureCanvasWxAgg(panel1, -1, self.figure1c)
        hbox1.Add(self.canvas1c, 1, wx.EXPAND)
        panel1.SetSizer(hbox1)
        vbox.Add(panel1, 1, wx.EXPAND)
        
        # Add a time vs DM vs SNR (with toolbar)
        panel3 = wx.Panel(self, -1)
        hbox3 = wx.BoxSizer(wx.VERTICAL)
        self.figure2 = Figure(figsize=(6,2))
        self.canvas2 = FigureCanvasWxAgg(panel3, -1, self.figure2)
        self.toolbar = RefreshAwareToolbar(self.canvas2, refreshCallback=lambda: self.data.draw(is_callback=True))
        self.toolbar.Realize()
        hbox3.Add(self.canvas2, 1, wx.ALIGN_LEFT | wx.EXPAND)
        hbox3.Add(self.toolbar, 0, wx.ALIGN_LEFT)
        panel3.SetSizer(hbox3)
        vbox.Add(panel3, 1, wx.EXPAND)
        self.panel3 = panel3
        
        # Use some sizers to see layout options
        self.SetSizer(vbox)
        self.SetAutoLayout(1)
        vbox.Fit(self)
        
    def initEvents(self):
        self.Bind(wx.EVT_MENU, self.onOpen, id=ID_OPEN)
        self.Bind(wx.EVT_MENU, self.onExit, id=ID_QUIT)
        
        self.Bind(wx.EVT_MENU, self.onColorValue, id=ID_COLOR_VALUE_DM)
        self.Bind(wx.EVT_MENU, self.onColorValue, id=ID_COLOR_VALUE_SNR)
        self.Bind(wx.EVT_MENU, self.onColorValue, id=ID_COLOR_VALUE_TIME)
        self.Bind(wx.EVT_MENU, self.onColorValue, id=ID_COLOR_VALUE_WIDTH)
        self.Bind(wx.EVT_MENU, self.onColorMap, id=ID_COLOR_MAP_PAIRED)
        self.Bind(wx.EVT_MENU, self.onColorMap, id=ID_COLOR_MAP_SPECTRAL)
        self.Bind(wx.EVT_MENU, self.onColorMap, id=ID_COLOR_MAP_BONE)
        self.Bind(wx.EVT_MENU, self.onColorMap, id=ID_COLOR_MAP_JET)
        self.Bind(wx.EVT_MENU, self.onColorMap, id=ID_COLOR_MAP_EARTH)
        self.Bind(wx.EVT_MENU, self.onColorMap, id=ID_COLOR_MAP_HEAT)
        self.Bind(wx.EVT_MENU, self.onColorMap, id=ID_COLOR_MAP_NCAR)
        self.Bind(wx.EVT_MENU, self.onColorMap, id=ID_COLOR_MAP_RAINBOW)
        self.Bind(wx.EVT_MENU, self.onColorMap, id=ID_COLOR_MAP_STERN)
        self.Bind(wx.EVT_MENU, self.onColorMap, id=ID_COLOR_MAP_GRAY)
        self.Bind(wx.EVT_MENU, self.onColorMap, id=ID_COLOR_INVERT)
        self.Bind(wx.EVT_MENU, self.onColorStretch, id=ID_COLOR_STRETCH_LINEAR)
        self.Bind(wx.EVT_MENU, self.onColorStretch, id=ID_COLOR_STRETCH_LOG)
        self.Bind(wx.EVT_MENU, self.onColorStretch, id=ID_COLOR_STRETCH_SQRT)
        self.Bind(wx.EVT_MENU, self.onColorStretch, id=ID_COLOR_STRETCH_SQRD)
        self.Bind(wx.EVT_MENU, self.onColorStretch, id=ID_COLOR_STRETCH_ASINH)
        self.Bind(wx.EVT_MENU, self.onColorStretch, id=ID_COLOR_STRETCH_SINH)
        self.Bind(wx.EVT_MENU, self.onColorStretch, id=ID_COLOR_STRETCH_HIST)
        
        self.Bind(wx.EVT_MENU, self.onDataSymbol, id=ID_DATA_SYMBOL_CIRCLE)
        self.Bind(wx.EVT_MENU, self.onDataSymbol, id=ID_DATA_SYMBOL_SQUARE)
        self.Bind(wx.EVT_MENU, self.onDataSymbol, id=ID_DATA_SYMBOL_DIAMOND)
        self.Bind(wx.EVT_MENU, self.onDataSymbol, id=ID_DATA_SYMBOL_HEXAGON)
        self.Bind(wx.EVT_MENU, self.onDataSymbol, id=ID_DATA_SYMBOL_PLUS)
        self.Bind(wx.EVT_MENU, self.onDataSize, id=ID_DATA_SIZE_DM)
        self.Bind(wx.EVT_MENU, self.onDataSize, id=ID_DATA_SIZE_SNR)
        self.Bind(wx.EVT_MENU, self.onDataSize, id=ID_DATA_SIZE_TIME)
        self.Bind(wx.EVT_MENU, self.onDataSize, id=ID_DATA_SIZE_WIDTH)
        self.Bind(wx.EVT_MENU, self.onDataAdjust, id=ID_DATA_ADJUST)
        
        self.Bind(wx.EVT_MENU, self.onDisplayDecimate, id=ID_DISPLAY_DECIMATE)
        self.Bind(wx.EVT_MENU, self.onDisplayDecimateAdjust, id=ID_DISPLAY_DECIMATE_ADJUST)
        
        self.Bind(wx.EVT_MENU, self.onHelp, id=ID_HELP)
        self.Bind(wx.EVT_MENU, self.onAbout, id=ID_ABOUT)
        
        # Key events
        self.canvas1a.Bind(wx.EVT_KEY_UP, self.onKeyPress)
        self.canvas1b.Bind(wx.EVT_KEY_UP, self.onKeyPress)
        self.canvas2.Bind(wx.EVT_KEY_UP,  self.onKeyPress)
        
        # Make the plots re-sizable
        # self.Bind(wx.EVT_PAINT, self.resizePlots)
        self.Bind(wx.EVT_SIZE, self.onSize)
        
        # Window manager close
        self.Bind(wx.EVT_CLOSE, self.onExit)
        
    def onOpen(self, event):
        """
        Open a file.
        """
        
        dlg = wx.FileDialog(self, "Choose a file", self.dirname, "", "Single Pulse (*.singlePulse)|*.singlePulse|All Files|*.*", wx.FD_OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            self.filename = dlg.GetFilename()
            self.dirname = dlg.GetDirectory()
            self.data = SinglePulse_GUI(self)
            
        dlg.Destroy()
        
    def onHelp(self, event):
        """
        Display the help window.
        """
        
        HelpWindow(self)
        
    def onAbout(self, event):
        """
        Display a very very very brief 'about' window.
        """
        
        dialog = wx.AboutDialogInfo()
        
        #dialog.SetIcon(wx.Icon(os.path.join(self.scriptPath, 'icons', 'lwa.png'), wx.BITMAP_TYPE_PNG))
        dialog.SetName('plotSinglePulse')
        dialog.SetVersion(__version__)
        dialog.SetDescription("""GUI for plotting single pulse data from PRESTO.\n\nLSL Version: %s""" % lsl.version.version)
        dialog.SetWebSite('http://lwa.unm.edu')
        dialog.AddDeveloper(__author__)
        
        # Debuggers/testers
        dialog.AddDocWriter(__author__)
        dialog.AddDocWriter("Kevin Stovall")
        
        wx.AboutBox(dialog)
        
    def onExit(self, event):
        """
        Quit plotWaterfall.
        """
        
        self.Destroy()
        
    def onColorValue(self, event):
        """
        Set the color coding to the specified quantity and refresh the plots.
        """
        
        wx.BeginBusyCursor()
        wx.Yield()
        
        if self.vmapMenu.IsChecked(ID_COLOR_VALUE_DM):
            idx = 0
        elif self.vmapMenu.IsChecked(ID_COLOR_VALUE_SNR):
            idx = 1
        elif self.vmapMenu.IsChecked(ID_COLOR_VALUE_TIME):
            idx = 2
        else:
            idx = 4
            
        if self.data.colorProperty != idx:
            self.data.colorProperty = idx
            self.data.draw()
            
        wx.EndBusyCursor()
        
    def onColorMap(self, event):
        """
        Set the colormap to the specified value and refresh the plots.
        """
        
        wx.BeginBusyCursor()
        wx.Yield()
        
        if self.cmapMenu.IsChecked(ID_COLOR_MAP_PAIRED):
            name = 'Paired'
        elif self.cmapMenu.IsChecked(ID_COLOR_MAP_SPECTRAL):
            name = 'Spectral'
        elif self.cmapMenu.IsChecked(ID_COLOR_MAP_BONE):
            name = 'bone'
        elif self.cmapMenu.IsChecked(ID_COLOR_MAP_EARTH):
            name = 'gist_earth'
        elif self.cmapMenu.IsChecked(ID_COLOR_MAP_HEAT):
            name = 'gist_heat'
        elif self.cmapMenu.IsChecked(ID_COLOR_MAP_NCAR):
            name = 'gist_ncar'
        elif self.cmapMenu.IsChecked(ID_COLOR_MAP_RAINBOW):
            name = 'gist_rainbow'
        elif self.cmapMenu.IsChecked(ID_COLOR_MAP_STERN):
            name = 'gist_stern'
        elif self.cmapMenu.IsChecked(ID_COLOR_MAP_GRAY):
            name = 'gist_gray'
        else:
            name = 'jet'
            
        # Check for inversion
        if self.cmapMenu.IsChecked(ID_COLOR_INVERT):
            if name == 'gist_gray':
                name = 'gist_yarg'
            else:
                name = '%s_r' % name
                
        newCM = cm.get_cmap(name)
        if self.data.cmap != newCM:
            self.data.cmap = newCM
            self.data.draw()
            
        wx.EndBusyCursor()
        
    def onColorStretch(self, event):
        """
        Set the color stretch to the specified scheme and refresh the plots.
        """
        
        wx.BeginBusyCursor()
        wx.Yield()
        
        if self.smapMenu.IsChecked(ID_COLOR_STRETCH_LOG):
            newN = LogNorm
        elif self.smapMenu.IsChecked(ID_COLOR_STRETCH_SQRT):
            newN = SqrtNorm
        elif self.smapMenu.IsChecked(ID_COLOR_STRETCH_SQRD):
            newN = SqrdNorm
        elif self.smapMenu.IsChecked(ID_COLOR_STRETCH_ASINH):
            newN = AsinhNorm
        elif self.smapMenu.IsChecked(ID_COLOR_STRETCH_SINH):
            newN = SinhNorm
        elif self.smapMenu.IsChecked(ID_COLOR_STRETCH_HIST):
            newN = HistEqNorm
        else:
            newN = Normalize
            
        if self.data.norm != newN:
            self.data.norm = newN
            self.data.draw()
            
        wx.EndBusyCursor()
        
    def onDataSymbol(self, event):
        """
        Set the plotting symbol to the specified shape and refresh the plots.
        """
        
        wx.BeginBusyCursor()
        wx.Yield()
        
        if self.mmapMenu.IsChecked(ID_DATA_SYMBOL_SQUARE):
            marker = 's'
        elif self.mmapMenu.IsChecked(ID_DATA_SYMBOL_DIAMOND):
            marker = 'D'
        elif self.mmapMenu.IsChecked(ID_DATA_SYMBOL_HEXAGON):
            marker = 'h'
        elif self.mmapMenu.IsChecked(ID_DATA_SYMBOL_PLUS):
            marker = '+'
        else:
            marker = 'o'
            
        if self.data.plotSymbol != marker:
            self.data.plotSymbol = marker
            self.data.draw()
            
        wx.EndBusyCursor()
        
    def onDataSize(self, event):
        """
        Set the symbol size coding to the specified quantity and refresh the plots.
        """
        
        wx.BeginBusyCursor()
        wx.Yield()
        
        if self.amapMenu.IsChecked(ID_DATA_SIZE_DM):
            idx = 0
        elif self.amapMenu.IsChecked(ID_DATA_SIZE_SNR):
            idx = 1
        elif self.amapMenu.IsChecked(ID_DATA_SIZE_TIME):
            idx = 2
        else:
            idx = 4
            
        if self.data.sizeProperty != idx:
            self.data.sizeProperty = idx
            self.data.draw()
            
        wx.EndBusyCursor()
        
    def onDataAdjust(self, event):
        """
        Bring up the data threshold adjustment dialog window.
        """
        
        ThresholdAdjust(self)
        
    def onDisplayDecimate(self, event):
        """
        Toggle scatter plot display decimation on/off.
        """
        
        wx.BeginBusyCursor()
        wx.Yield()
        
        self.data.fullRes = not self.data.fullRes
        self.dadj.Enable(not self.data.fullRes)
        self.data.draw()
        
        wx.EndBusyCursor()
        
    def onDisplayDecimateAdjust(self, event):
        """
        Bring up the scatter plot display decimation adjustment dialog window.
        """
        
        DecimationAdjust(self)
        
    def onKeyPress(self, event):
        """
        Move the current spectra mark up and down with a keypress.
        """
        
        keycode = event.GetKeyCode()
        if keycode < 256:
            keycode = chr(keycode)
            
        if keycode == wx.WXK_UP or keycode == wx.WXK_RIGHT:
            ## Move forward one pulse
            if self.data.pulseClick is not None:
                pass
            
        elif keycode == wx.WXK_DOWN or keycode == wx.WXK_LEFT:
            ## Move backward one pulse
            if self.data.pulseClick is not None:
                pass
            
        elif keycode == 'M':
            pass
            
        elif keycode == 'U':
            pass
            
        else:
            pass
            
        event.Skip()
        
    def onSize(self, event):
        event.Skip()
        
        wx.CallAfter(self.resizePlots)
        
    def resizePlots(self, event=None):
        # Get the current size of the window and the navigation toolbar
        w, h = self.GetClientSize()
        wt, ht = self.toolbar.GetSize()
        
        # Come up with new figure sizes in inches.
        # NOTE:  The height of the waterfall plot at the bottom needs to be
        #        corrected for the height of the toolbar
        dpi = self.figure1a.get_dpi()
        newW0 = 1.0*w/dpi
        newW1 = 1.0*(2.*w/6)/dpi
        newW2 = 1.0*(2.*w/6)/dpi
        newW3 = 1.0*(2.*w/6)/dpi
        newH0 = 1.0*(h/2)/dpi
        newH1 = 1.0*(h/2-ht)/dpi
        
        # Set the figure sizes and redraw
        self.figure1a.set_size_inches((newW1, newH0))
        self.figure1a.tight_layout()
        self.figure1a.canvas.draw()
        
        self.figure1b.set_size_inches((newW2, newH0))
        self.figure1b.tight_layout()
        self.figure1b.canvas.draw()
        
        self.figure1c.set_size_inches((newW3, newH0))
        self.figure1c.tight_layout()
        self.figure1c.canvas.draw()
        
        self.figure2.set_size_inches((newW0, newH1))
        self.figure2.tight_layout()
        self.figure2.canvas.draw()
        
        self.panel3.Refresh()
        
    def GetToolBar(self):
        # You will need to override GetToolBar if you are using an 
        # unmanaged toolbar in your frame
        return self.toolbar


ID_THRESHOLD_THRESH_INC = 100
ID_THRESHOLD_THRESH_DEC = 101
ID_THRESHOLD_WIDTH_MIN_INC = 102
ID_THRESHOLD_WIDTH_MIN_DEC = 103
ID_THRESHOLD_WIDTH_MAX_INC = 104
ID_THRESHOLD_WIDTH_MAX_DEC = 105
ID_THRESHOLD_OK = 106

class ThresholdAdjust(wx.Frame):
    def __init__ (self, parent):	
        wx.Frame.__init__(self, parent, title='Pulse Filtering Adjustment')
        
        self.parent = parent
        
        self.initUI()
        self.initEvents()
        self.Show()
        
    def initUI(self):
        row = 0
        panel = wx.Panel(self)
        sizer = wx.GridBagSizer(0, 0)
        
        thr = wx.StaticText(panel, label='Min. S/N Threshold:')
        thrText = wx.TextCtrl(panel)
        thrDec = wx.Button(panel, ID_THRESHOLD_THRESH_DEC, '-', size=(56, 28))
        thrInc = wx.Button(panel, ID_THRESHOLD_THRESH_INC, '+', size=(56, 28))
        
        lwr = wx.StaticText(panel, label='Min. Pulse Width [ms]:')
        lwrText = wx.TextCtrl(panel)
        lwrDec = wx.Button(panel, ID_THRESHOLD_WIDTH_MIN_DEC, '-', size=(56, 28))
        lwrInc = wx.Button(panel, ID_THRESHOLD_WIDTH_MIN_INC, '+', size=(56, 28))
        
        upr = wx.StaticText(panel, label='Max. Pulse Width [ms]:')
        uprText = wx.TextCtrl(panel)
        uprDec = wx.Button(panel, ID_THRESHOLD_WIDTH_MAX_DEC, '-', size=(56, 28))
        uprInc = wx.Button(panel, ID_THRESHOLD_WIDTH_MAX_INC, '+', size=(56, 28))
        
        sizer.Add(thr,     pos=(row+0, 0), span=(1, 1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
        sizer.Add(thrText, pos=(row+0, 1), span=(1, 1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
        sizer.Add(thrDec,  pos=(row+0, 2), span=(1, 1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
        sizer.Add(thrInc,  pos=(row+0, 3), span=(1, 1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
        sizer.Add(upr,     pos=(row+1, 0), span=(1, 1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
        sizer.Add(uprText, pos=(row+1, 1), span=(1, 1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
        sizer.Add(uprDec,  pos=(row+1, 2), span=(1, 1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
        sizer.Add(uprInc,  pos=(row+1, 3), span=(1, 1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
        sizer.Add(lwr,     pos=(row+2, 0), span=(1, 1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
        sizer.Add(lwrText, pos=(row+2, 1), span=(1, 1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
        sizer.Add(lwrDec,  pos=(row+2, 2), span=(1, 1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
        sizer.Add(lwrInc,  pos=(row+2, 3), span=(1, 1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
        
        line = wx.StaticLine(panel)
        sizer.Add(line, pos=(row+3, 0), span=(1, 4), flag=wx.EXPAND|wx.BOTTOM, border=10)
        
        row += 4
        
        #
        # Buttons
        #
        
        ok = wx.Button(panel, ID_THRESHOLD_OK, 'Ok', size=(56, 28))
        sizer.Add(ok, pos=(row+0, 3), flag=wx.ALL, border=5)
        
        panel.SetSizer(sizer)
        sizer.Fit(self)

        self.tText = thrText
        self.lText = lwrText
        self.uText = uprText
        
        #
        # Set current values
        #
        self.tText.SetValue('%.1f' % self.parent.data.dataThreshold[0])
        self.lText.SetValue('%.3f' % self.parent.data.dataThreshold[1])
        self.uText.SetValue('%.3f' % self.parent.data.dataThreshold[2])
        
    def initEvents(self):
        self.Bind(wx.EVT_BUTTON, self.onThresholdDecrease, id=ID_THRESHOLD_THRESH_DEC)
        self.Bind(wx.EVT_BUTTON, self.onThresholdIncrease, id=ID_THRESHOLD_THRESH_INC)
        self.Bind(wx.EVT_BUTTON, self.onUpperDecrease, id=ID_THRESHOLD_WIDTH_MAX_DEC)
        self.Bind(wx.EVT_BUTTON, self.onUpperIncrease, id=ID_THRESHOLD_WIDTH_MAX_INC)
        self.Bind(wx.EVT_BUTTON, self.onLowerDecrease, id=ID_THRESHOLD_WIDTH_MIN_DEC)
        self.Bind(wx.EVT_BUTTON, self.onLowerIncrease, id=ID_THRESHOLD_WIDTH_MIN_INC)
        
        self.Bind(wx.EVT_BUTTON, self.onOk, id=ID_THRESHOLD_OK)
        
        self.Bind(wx.EVT_KEY_UP, self.onKeyPress)
        
    def onKeyPress(self, event):
        keycode = event.GetKeyCode()
        
        if keycode == wx.WXK_RETURN:
            self.parent.data.dataThreshold[0] = float(self.tText.GetValue())
            self.parent.data.dataThreshold[1] = float(self.lText.GetValue())
            self.parent.data.dataThreshold[2] = float(self.uText.GetValue())
            self.parent.data.draw(recompute=True)
        else:
            pass
        
    def onThresholdDecrease(self, event):
        self.parent.data.dataThreshold[0] -= 1
        if self.parent.data.dataThreshold[0] < 1.0:
            self.parent.data.dataThreshold[0] = 1.0
        self.tText.SetValue('%.1f' % self.parent.data.dataThreshold[0])
        self.parent.data.draw(recompute=True)
        
    def onThresholdIncrease(self, event):
        self.parent.data.dataThreshold[0] += 1
        self.tText.SetValue('%.1f' % self.parent.data.dataThreshold[0])
        self.parent.data.draw(recompute=True)
        
    def onUpperDecrease(self, event):
        self.parent.data.dataThreshold[2] -= self.parent.data.dataThreshold[2]*0.1
        if self.parent.data.dataThreshold[2] < 0.0:
            self.parent.data.dataThreshold[2] = 0.0
        self.uText.SetValue('%.3f' % self.parent.data.dataThreshold[2])
        self.parent.data.draw(recompute=True)
        
    def onUpperIncrease(self, event):
        self.parent.data.dataThreshold[2] += self.parent.data.dataThreshold[2]*0.1
        self.uText.SetValue('%.3f' % self.parent.data.dataThreshold[2])
        self.parent.data.draw(recompute=True)
        
    def onLowerDecrease(self, event):
        self.parent.data.dataThreshold[1] -= self.parent.data.dataThreshold[1]*0.1
        if self.parent.data.dataThreshold[1] < 0.0:
            self.parent.data.dataThreshold[1] = 0.0
        self.lText.SetValue('%.3f' % self.parent.data.dataThreshold[1])
        self.parent.data.draw(recompute=True)
        
    def onLowerIncrease(self, event):
        self.parent.data.dataThreshold[1] += self.parent.data.dataThreshold[1]*0.1
        self.lText.SetValue('%.3f' % self.parent.data.dataThreshold[1])
        self.parent.data.draw(recompute=True)
        
    def onOk(self, event):
        needToRedraw = False
        
        if self.parent.data.dataThreshold[0] != float(self.tText.GetValue()):
            self.parent.data.dataThreshold[0] = float(self.tText.GetValue())
            needToRedraw = True
        if self.parent.data.dataThreshold[1] != float(self.lText.GetValue()):
            self.parent.data.dataThreshold[1] = float(self.lText.GetValue())
            needToRedraw = True
        if self.parent.data.dataThreshold[2] != float(self.uText.GetValue()):
            self.parent.data.dataThreshold[2] = float(self.uText.GetValue())
            needToRedraw = True
            
        if needToRedraw:
            self.parent.data.draw(recompute=True)
            
        self.Close()


ID_DECIMATION_INC = 200
ID_DECIMATION_DEC = 201
ID_DECIMATION_OK = 202

class DecimationAdjust(wx.Frame):
    def __init__ (self, parent):	
        wx.Frame.__init__(self, parent, title='Display Decimation Adjustment')
        
        self.parent = parent
        
        self.initUI()
        self.initEvents()
        self.Show()
        
    def initUI(self):
        row = 0
        panel = wx.Panel(self)
        sizer = wx.GridBagSizer(0, 0)
        
        upr = wx.StaticText(panel, label='Max. Points to Plot:')
        uprText = wx.TextCtrl(panel)
        uprDec = wx.Button(panel, ID_DECIMATION_DEC, '-', size=(56, 28))
        uprInc = wx.Button(panel, ID_DECIMATION_INC, '+', size=(56, 28))
        
        sizer.Add(upr,     pos=(row+0, 0), span=(1, 1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
        sizer.Add(uprText, pos=(row+0, 1), span=(1, 1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
        sizer.Add(uprDec,  pos=(row+0, 2), span=(1, 1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
        sizer.Add(uprInc,  pos=(row+0, 3), span=(1, 1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
        
        line = wx.StaticLine(panel)
        sizer.Add(line, pos=(row+1, 0), span=(1, 4), flag=wx.EXPAND|wx.BOTTOM, border=10)
        
        row += 2
        
        #
        # Buttons
        #
        
        ok = wx.Button(panel, ID_DECIMATION_OK, 'Ok', size=(56, 28))
        sizer.Add(ok, pos=(row+0, 3), flag=wx.ALL, border=5)
        
        panel.SetSizer(sizer)
        sizer.Fit(self)

        self.uText = uprText
        
        #
        # Set current values
        #
        self.uText.SetValue('%i' % self.parent.data.maxPoints)
        
    def initEvents(self):
        self.Bind(wx.EVT_BUTTON, self.onUpperDecrease, id=ID_DECIMATION_DEC)
        self.Bind(wx.EVT_BUTTON, self.onUpperIncrease, id=ID_DECIMATION_INC)
        
        self.Bind(wx.EVT_BUTTON, self.onOk, id=ID_DECIMATION_OK)
        
        self.Bind(wx.EVT_KEY_UP, self.onKeyPress)
        
    def onKeyPress(self, event):
        keycode = event.GetKeyCode()
        
        if keycode == wx.WXK_RETURN:
            self.parent.data.maxPoints = int(self.uText.GetValue())
            self.parent.data.draw(recompute=True)
        else:
            pass
        
    def onUpperDecrease(self, event):
        self.parent.data.maxPoints -= 1000
        if self.parent.data.maxPoints < 1000:
            self.parent.data.maxPoints = 1000
        self.uText.SetValue('%i' % self.parent.data.maxPoints)
        
        t0 = time.time()
        self.parent.data.draw(recompute=True)
        t1 = time.time()
        print("-> Drawing time with %i points is %.3f s" % (self.parent.data.maxPoints, t1-t0))
        
    def onUpperIncrease(self, event):
        self.parent.data.maxPoints += 1000
        self.uText.SetValue('%i' % self.parent.data.maxPoints)
        
        t0 = time.time()
        self.parent.data.draw(recompute=True)
        t1 = time.time()
        print("-> Drawing time with %i points is %.3f s" % (self.parent.data.maxPoints, t1-t0))
        
    def onOk(self, event):
        needToRedraw = False
        
        if self.parent.data.maxPoints != int(self.uText.GetValue()):
            self.parent.data.maxPoints = int(self.uText.GetValue())
            needToRedraw = True
            
        if needToRedraw:
            self.parent.data.draw(recompute=True)
            
        self.Close()


class SliceDisplay(wx.Frame):
    """
    Window for displaying a time slice of data in a zoomable fashion
    """
    
    def __init__(self, parent, t, dm, width):
        wx.Frame.__init__(self, parent, title='Time Slice', size=(400, 375))
        
        self.parent = parent
        
        self.t = t
        self.dm = dm
        self.width = width
        
        self.initUI()
        self.initEvents()
        self.Show()
        
        self.initPlot()
        
    def initUI(self):
        """
        Start the user interface.
        """
        
        self.statusbar = self.CreateStatusBar()
        
        hbox = wx.BoxSizer(wx.HORIZONTAL)
        
        # Add plots to panel 1
        panel1 = wx.Panel(self, -1)
        vbox1 = wx.BoxSizer(wx.VERTICAL)
        self.figure = Figure()
        self.canvas = FigureCanvasWxAgg(panel1, -1, self.figure)
        self.toolbar = NavigationToolbar2WxAgg(self.canvas)
        self.toolbar.Realize()
        vbox1.Add(self.canvas,  1, wx.ALIGN_LEFT | wx.EXPAND)
        vbox1.Add(self.toolbar, 0, wx.ALIGN_LEFT)
        panel1.SetSizer(vbox1)
        hbox.Add(panel1, 1, wx.EXPAND)
        
        # Use some sizers to see layout options
        self.SetSizer(hbox)
        self.SetAutoLayout(1)
        hbox.Fit(self)
        
    def initEvents(self):
        """
        Set all of the various events in the data range window.
        """
        
        # Make the images resizable
        self.Bind(wx.EVT_PAINT, self.resizePlots)
        
    def initPlot(self):
        """
        Populate the figure/canvas areas with a plot.  We only need to do this
        once for this type of window.
        """
        
        # Plot the slice
        self.figure.clf()
        self.ax1 = self.figure.gca()
        
        ## Determine what data to plot
        fLow, fHigh = self.parent.data.meta.lofreq, self.parent.data.meta.lofreq + self.parent.data.meta.BW
        
        slope = -_D*(1.0/fLow**2 - 1.0/fHigh**2)
        
        sLow, wLow, wHigh = self.parent.data.dataThreshold
        tCutLow = self.t + (self.parent.data.dmMax-self.dm)*slope - 1
        tCutHigh = self.t + (self.parent.data.dmMin-self.dm)*slope + 1
        
        valid = numpy.where( (self.parent.data.data[:,2] >= tCutLow) & \
                        (self.parent.data.data[:,2] <= tCutHigh) & \
                        (self.parent.data.data[:,1] >= sLow ) & \
                        (self.parent.data.data[:,4] >= wLow ) & \
                        (self.parent.data.data[:,4] <= wHigh ) )[0]
        if len(valid) == 0:
            self.ax1.set_title('Nothing to display')
            return False
            
        self.subData = self.parent.data.data[valid,:]
        
        deltaT = self.subData[:,2] + (self.subData[:,0]-self.dm)*slope - self.t
        valid = numpy.where( numpy.abs(deltaT) < 1 )[0]
        if len(valid) == 0:
            self.ax1.set_title('Nothing to display')
            return False
            
        self.subData = self.subData[valid,:]
        
        m = self.ax1.scatter(self.subData[:,0], self.subData[:,1], 
                            c=self.subData[:,self.parent.data.colorProperty], 
                            s=self.subData[:,self.parent.data.sizeProperty]*5, 
                            cmap=self.parent.data.cmap, norm=self.parent.data.norm(*self.parent.data.limits[self.parent.data.colorProperty]), 
                            marker='o', edgecolors='face')
        try:
            cm = self.frame.figure.colorbar(m, use_gridspec=True)
        except:
            if len(self.figure.get_axes()) > 1:
                self.figure.delaxes( self.figure.get_axes()[-1] )
            cm = self.figure.colorbar(m)
        if self.parent.data.colorProperty == 0:
            cm.ax.set_ylabel('DM [pc cm$^{-3}$]')
        elif self.parent.data.colorProperty == 1:
            cm.ax.set_ylabel('S/N')
        elif self.parent.data.colorProperty == 2:
            cm.ax.set_ylabel('Elapsed Time [s]')
        else:
            cm.ax.set_ylabel('Width [ms]')
            
        ## Generate the "expected S/N vs. DM curve from Cordes & McLaughlin (2003, Apj, 596, 1142)
        best = numpy.argmax(self.subData[:,1])
        dm = self.subData[best,0]
        width = self.subData[best,4]
        
        dms1 = numpy.linspace(self.subData[:,0].min(), self.subData[:,0].max(), num=101)
        dms0 = numpy.linspace(self.parent.data.ax2.get_ylim()[0], dms1.min(), num=101)
        dms2 = numpy.linspace(dms1.max(), self.parent.data.ax2.get_ylim()[1], num=101)
        dms = numpy.concatenate([dms0, dms1, dms2])
        zeta = 6.91e-3 * (dms - dm) * self.parent.data.meta.BW / width / ((fLow+fHigh)/2.0/1000.0)**3
        thr = numpy.sqrt(numpy.pi)/2 * erf(zeta)/zeta
        thr *= self.subData[:,1].max()
        self.ax1.plot(dms, thr, linestyle='--', marker='', color='black', label='Cordes & McLaughlin (2003)')
        self.ax1.legend(loc=0)
        
        self.ax1.axis('auto')
        self.ax1.set_xlim(self.parent.data.ax2.get_ylim())
        self.ax1.set_ylim((self.subData[:,1].min()*0.90, self.subData[:,1].max()*1.1))
        self.ax1.set_xlabel('DM [pc cm$^{-3}$]')
        self.ax1.set_ylabel('S/N')
        self.ax1.set_title('Slice through %.1f s, %.3f pc cm$^{-3}$' % (self.t, self.dm))
        
        ## Draw and save the click (Why?)
        self.canvas.draw()
        self.connect()
        
    def connect(self):
        """
        Connect to all the events we need to interact with the plots.
        """
        
        self.cidmotion  = self.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)
    
    def on_motion(self, event):
        """
        Deal with motion events in the stand field window.  This involves 
        setting the status bar with the current x and y coordinates as well
        as the stand number of the selected stand (if any).
        """
        
        if event.inaxes:
            clickX = event.xdata
            clickY = event.ydata
            
            d = (self.subData[:,0] - clickX)**2 + (self.subData[:,1] - clickY)**2
            best = numpy.where( (d == d.min()) & \
                            (self.subData[:,4] >= self.parent.data.dataThreshold[1]) & \
                            (self.subData[:,4] <= self.parent.data.dataThreshold[2]) )[0][0]
            
            dm, snr, width = self.subData[best,0], self.subData[best,1], self.subData[best,4]
            
            self.statusbar.SetStatusText("DM=%.4f pc cm^-3, S/N=%.1f, width=%.3f ms" % (dm, snr, width))
        else:
            self.statusbar.SetStatusText("")
    
    def disconnect(self):
        """
        Disconnect all the stored connection ids.
        """
        
        self.figure.canvas.mpl_disconnect(self.cidmotion)
        
    def onCancel(self, event):
        self.Close()
        
    def resizePlots(self, event):
        # Get the current size of the window and the navigation toolbar
        w, h = self.GetClientSize()
        wt, ht = self.toolbar.GetSize()
        
        dpi = self.figure.get_dpi()
        newW = 1.0*w/dpi
        newH = 1.0*(h-ht)/dpi
        self.figure.set_size_inches((newW, newH))
        self.figure.tight_layout()
        self.figure.canvas.draw()

    def GetToolBar(self):
        # You will need to override GetToolBar if you are using an 
        # unmanaged toolbar in your frame
        return self.toolbar


ID_WATERFALL_AUTO = 30000
ID_WATERFALL_ADJUST = 30001
ID_WATERFALL_MAP_PAIRED = 30002
ID_WATERFALL_MAP_SPECTRAL = 30003
ID_WATERFALL_MAP_BONE = 30004
ID_WATERFALL_MAP_JET = 30005
ID_WATERFALL_MAP_EARTH = 30006
ID_WATERFALL_MAP_HEAT = 30007
ID_WATERFALL_MAP_NCAR = 30008
ID_WATERFALL_MAP_RAINBOW = 30009
ID_WATERFALL_MAP_STERN = 30010
ID_WATERFALL_MAP_GRAY = 30011
ID_WATERFALL_INVERT = 30012
ID_WATERFALL_STRETCH_LINEAR = 30013
ID_WATERFALL_STRETCH_LOG = 30014
ID_WATERFALL_STRETCH_SQRT = 30015
ID_WATERFALL_STRETCH_SQRD = 30016
ID_WATERFALL_STRETCH_ASINH = 30017
ID_WATERFALL_STRETCH_SINH = 30018
ID_WATERFALL_STRETCH_HIST = 30019

ID_WATERFALL_PRODUCT_1 = 30020
ID_WATERFALL_DECIMATION_AUTO = 30030
ID_WATERFALL_DECIMATION_ADJUST = 30031

ID_WATERFALL_BANDPASS_ON = 30032
ID_WATERFALL_BANDPASS_OFF = 30033

ID_WATERFALL_SWEEP = 30034
ID_WATERFALL_PROFILE = 30035


class WaterfallDisplay(wx.Frame):
    """
    Window for displaying the waterfall data in a zoomable fashion
    """
    
    def __init__(self, parent, fitsname, t, dm, width):
        wx.Frame.__init__(self, parent, title='Waterfall', size=(400, 375))
        
        self.parent = parent
        self.fitsname = fitsname
        self.t = t
        self.dm = dm
        self.width = width / 1000.0	# ms -> s
        
        # Convert time from barycentric to topocentric, if required
        if self.parent.data.meta.bary:
            if self.parent.data.bary2topo is not None:
                self.t = self.parent.data.bary2topo(self.t)
                
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
        self.initEvents()
        self.Show()
        
        self.render()
        self.draw()
        
    def initUI(self):
        """
        Start the user interface.
        """
        
        self.statusbar = self.CreateStatusBar()
        
        menuBar = wx.MenuBar()
        colorMenu = wx.Menu()
        dataMenu = wx.Menu()
        bandpassMenu = wx.Menu()
        displayMenu = wx.Menu()
        
        ## Color Menu
        auto = wx.MenuItem(colorMenu, ID_WATERFALL_AUTO, '&Auto-scale Colorbar')
        AppendMenuItem(colorMenu, auto)
        cadj = wx.MenuItem(colorMenu, ID_WATERFALL_ADJUST, 'Adjust &Contrast')
        AppendMenuItem(colorMenu, cadj)
        cmap = wx.Menu()
        cmap.AppendRadioItem(ID_WATERFALL_MAP_PAIRED, '&Paired')
        cmap.AppendRadioItem(ID_WATERFALL_MAP_SPECTRAL, "&Spectral")
        cmap.AppendRadioItem(ID_WATERFALL_MAP_BONE, '&Bone')
        cmap.AppendRadioItem(ID_WATERFALL_MAP_JET, '&Jet')
        cmap.AppendRadioItem(ID_WATERFALL_MAP_EARTH, '&Earth')
        cmap.AppendRadioItem(ID_WATERFALL_MAP_HEAT, '&Heat')
        cmap.AppendRadioItem(ID_WATERFALL_MAP_NCAR, '&NCAR')
        cmap.AppendRadioItem(ID_WATERFALL_MAP_RAINBOW, '&Rainbow')
        cmap.AppendRadioItem(ID_WATERFALL_MAP_STERN, 'S&tern')
        cmap.AppendRadioItem(ID_WATERFALL_MAP_GRAY, '&Gray')
        cmap.AppendSeparator()
        cmap.AppendCheckItem(ID_WATERFALL_INVERT, 'In&vert')
        colorMenu.AppendMenu(-1, 'Color &Map', cmap)
        cmap.Check(ID_WATERFALL_MAP_JET, True)
        self.cmapMenu = cmap
        smap = wx.Menu()
        smap.AppendRadioItem(ID_WATERFALL_STRETCH_LINEAR, '&Linear')
        smap.AppendRadioItem(ID_WATERFALL_STRETCH_LOG, 'Lo&g')
        smap.AppendRadioItem(ID_WATERFALL_STRETCH_SQRT, 'Square &Root')
        smap.AppendRadioItem(ID_WATERFALL_STRETCH_SQRD, '&Squared')
        smap.AppendRadioItem(ID_WATERFALL_STRETCH_ASINH, '&ASinh')
        smap.AppendRadioItem(ID_WATERFALL_STRETCH_SINH, '&Sinh')
        smap.AppendRadioItem(ID_WATERFALL_STRETCH_HIST, '&Histogram Equalization')
        colorMenu.AppendMenu(-1, 'Color &Stretch', smap)
        smap.Check(ID_WATERFALL_STRETCH_LINEAR, True)
        self.smapMenu = smap
        
        ## Data Menu
        self.dataMenuOptions = []
        for i,dataProduct in enumerate(self.data_products):
            dmo = dataMenu.AppendRadioItem(ID_WATERFALL_PRODUCT_1+i, '%s' % dataProduct)
            self.dataMenuOptions.append( dmo )
        dataMenu.AppendSeparator()
        adec = dataMenu.AppendCheckItem(ID_WATERFALL_DECIMATION_AUTO, '&Auto Time Decimation')
        dataMenu.Check(ID_WATERFALL_DECIMATION_AUTO, True)
        self.adec = adec
        dadj = wx.MenuItem(dataMenu, ID_WATERFALL_DECIMATION_ADJUST, 'Adjust Time &Decimation')
        AppendMenuItem(dataMenu, dadj)
        dadj.Enable(False)
        self.dadj = dadj
            
        ## Bandpass Menu
        bandpassMenu.AppendRadioItem(ID_WATERFALL_BANDPASS_ON,  'On')
        bandpassMenu.AppendRadioItem(ID_WATERFALL_BANDPASS_OFF, 'Off')
        
        ## Display Menu
        sshw = displayMenu.AppendCheckItem(ID_WATERFALL_SWEEP, 'Show &Sweep')
        self.sshw = sshw
        pshw = displayMenu.AppendCheckItem(ID_WATERFALL_PROFILE, 'Show &Profile')
        self.pshw = pshw
        displayMenu.Check(ID_WATERFALL_SWEEP, True)
        displayMenu.Check(ID_WATERFALL_PROFILE, True)
        
        # Creating the menubar.
        menuBar.Append(colorMenu, "&Color")
        menuBar.Append(dataMenu, "&Data")
        menuBar.Append(bandpassMenu, "&Bandpass")
        menuBar.Append(displayMenu, "Di&splay")
        self.SetMenuBar(menuBar)  # Adding the MenuBar to the Frame content.
        
        # Save the various menus so that we can disable them if need be
        self.colorMenu = colorMenu
        self.dataMenu = dataMenu
        self.bandpassMenu = bandpassMenu
        self.displayMenu = displayMenu
        
        hbox = wx.BoxSizer(wx.HORIZONTAL)
        
        # Add plots to panel 1
        panel1 = wx.Panel(self, -1)
        vbox1 = wx.BoxSizer(wx.VERTICAL)
        self.figure = Figure()
        self.canvas = FigureCanvasWxAgg(panel1, -1, self.figure)
        self.toolbar = NavigationToolbar2WxAgg(self.canvas)
        self.toolbar.Realize()
        vbox1.Add(self.canvas,  1, wx.ALIGN_LEFT | wx.EXPAND)
        vbox1.Add(self.toolbar, 0, wx.ALIGN_LEFT)
        panel1.SetSizer(vbox1)
        hbox.Add(panel1, 1, wx.EXPAND)
        
        # Use some sizers to see layout options
        self.SetSizer(hbox)
        self.SetAutoLayout(1)
        hbox.Fit(self)
        
    def initEvents(self):
        """
        Set all of the various events in the data range window.
        """
        
        self.Bind(wx.EVT_MENU, self.onAutoscale, id=ID_WATERFALL_AUTO)
        self.Bind(wx.EVT_MENU, self.onAdjust, id=ID_WATERFALL_ADJUST)
        self.Bind(wx.EVT_MENU, self.onColorMap, id=ID_WATERFALL_MAP_PAIRED)
        self.Bind(wx.EVT_MENU, self.onColorMap, id=ID_WATERFALL_MAP_SPECTRAL)
        self.Bind(wx.EVT_MENU, self.onColorMap, id=ID_WATERFALL_MAP_BONE)
        self.Bind(wx.EVT_MENU, self.onColorMap, id=ID_WATERFALL_MAP_JET)
        self.Bind(wx.EVT_MENU, self.onColorMap, id=ID_WATERFALL_MAP_EARTH)
        self.Bind(wx.EVT_MENU, self.onColorMap, id=ID_WATERFALL_MAP_HEAT)
        self.Bind(wx.EVT_MENU, self.onColorMap, id=ID_WATERFALL_MAP_NCAR)
        self.Bind(wx.EVT_MENU, self.onColorMap, id=ID_WATERFALL_MAP_RAINBOW)
        self.Bind(wx.EVT_MENU, self.onColorMap, id=ID_WATERFALL_MAP_STERN)
        self.Bind(wx.EVT_MENU, self.onColorMap, id=ID_WATERFALL_MAP_GRAY)
        self.Bind(wx.EVT_MENU, self.onColorMap, id=ID_WATERFALL_INVERT)
        self.Bind(wx.EVT_MENU, self.onColorStretch, id=ID_WATERFALL_STRETCH_LINEAR)
        self.Bind(wx.EVT_MENU, self.onColorStretch, id=ID_WATERFALL_STRETCH_LOG)
        self.Bind(wx.EVT_MENU, self.onColorStretch, id=ID_WATERFALL_STRETCH_SQRT)
        self.Bind(wx.EVT_MENU, self.onColorStretch, id=ID_WATERFALL_STRETCH_SQRD)
        self.Bind(wx.EVT_MENU, self.onColorStretch, id=ID_WATERFALL_STRETCH_ASINH)
        self.Bind(wx.EVT_MENU, self.onColorStretch, id=ID_WATERFALL_STRETCH_SINH)
        self.Bind(wx.EVT_MENU, self.onColorStretch, id=ID_WATERFALL_STRETCH_HIST)
        
        for i,dataProduct in enumerate(self.data_products):
            self.Bind(wx.EVT_MENU, self.onDataProduct, id=ID_WATERFALL_PRODUCT_1+i)
        self.Bind(wx.EVT_MENU, self.onAutoDecimation, id=ID_WATERFALL_DECIMATION_AUTO)
        self.Bind(wx.EVT_MENU, self.onAdjustDecimation, id=ID_WATERFALL_DECIMATION_ADJUST)
            
        self.Bind(wx.EVT_MENU, self.onBandpassOn, id=ID_WATERFALL_BANDPASS_ON)
        self.Bind(wx.EVT_MENU, self.onBandpassOff, id=ID_WATERFALL_BANDPASS_OFF)
        
        self.Bind(wx.EVT_MENU, self.onShowSweep, id=ID_WATERFALL_SWEEP)
        self.Bind(wx.EVT_MENU, self.onShowProfile, id=ID_WATERFALL_PROFILE)
        
        # Make the images resizable
        self.Bind(wx.EVT_PAINT, self.resizePlots)
        
    def onAutoscale(self, event):
        """
        Auto-scale the current data display.
        """
        
        wx.BeginBusyCursor()
        wx.Yield()
        
        i = self.index
        toUse = numpy.arange(self.spec.shape[2]//10, 9*self.spec.shape[2]//10)
        
        try:
            from _helper import FastAxis1Percentiles5And99
            if self.bandpass:
                self.limitsBandpass[i] = list(FastAxis1Percentiles5And99(self.specBandpass.data, i, chanMin=self.spec.shape[2]//10, chanMax=9*self.spec.shape[2]//10))
            else:
                self.limits[i] = list(FastAxis1Percentiles5And99(self.spec.data, i))
        except ImportError:
            if self.bandpass:
                self.limitsBandpass[i] = numpy.percentile(self.specBandpass[:,i,toUse], (5, 99))
            else:
                self.limits[i] = numpy.percentile(self.spec[:,i,:], (5, 99))
                
        if self.usedB:
            if self.bandpass:
                self.limitsBandpass[i] = [to_dB(v) for v in self.limitsBandpass[i]]
            else:
                self.limits[i] = [to_dB(v) for v in self.limits[i]]
                
        self.draw()
        
        wx.EndBusyCursor()
        
    def onAdjust(self, event):
        """
        Bring up the colorbar adjustment dialog window.
        """
        
        WaterfallContrastAdjust(self)
        
    def onColorMap(self, event):
        """
        Set the colormap to the specified value and refresh the plots.
        """
        
        wx.BeginBusyCursor()
        wx.Yield()
        
        if self.cmapMenu.IsChecked(ID_WATERFALL_MAP_PAIRED):
            name = 'Paired'
        elif self.cmapMenu.IsChecked(ID_WATERFALL_MAP_SPECTRAL):
            name = 'Spectral'
        elif self.cmapMenu.IsChecked(ID_WATERFALL_MAP_BONE):
            name = 'bone'
        elif self.cmapMenu.IsChecked(ID_WATERFALL_MAP_EARTH):
            name = 'gist_earth'
        elif self.cmapMenu.IsChecked(ID_WATERFALL_MAP_HEAT):
            name = 'gist_heat'
        elif self.cmapMenu.IsChecked(ID_WATERFALL_MAP_NCAR):
            name = 'gist_ncar'
        elif self.cmapMenu.IsChecked(ID_WATERFALL_MAP_RAINBOW):
            name = 'gist_rainbow'
        elif self.cmapMenu.IsChecked(ID_WATERFALL_MAP_STERN):
            name = 'gist_stern'
        elif self.cmapMenu.IsChecked(ID_WATERFALL_MAP_GRAY):
            name = 'gist_gray'
        else:
            name = 'jet'
            
        # Check for inversion
        if self.cmapMenu.IsChecked(ID_WATERFALL_INVERT):
            if name == 'gist_gray':
                name = 'gist_yarg'
            else:
                name = '%s_r' % name
                
        newCM = cm.get_cmap(name)
        if self.cmap != newCM:
            self.cmap = newCM
            self.draw()
            
        wx.EndBusyCursor()
        
    def onColorStretch(self, event):
        """
        Set the color stretch to the specified scheme and refresh the plots.
        """
        
        wx.BeginBusyCursor()
        wx.Yield()
        
        if self.smapMenu.IsChecked(ID_WATERFALL_STRETCH_LOG):
            newN = LogNorm
        elif self.smapMenu.IsChecked(ID_WATERFALL_STRETCH_SQRT):
            newN = SqrtNorm
        elif self.smapMenu.IsChecked(ID_WATERFALL_STRETCH_SQRD):
            newN = SqrdNorm
        elif self.smapMenu.IsChecked(ID_WATERFALL_STRETCH_ASINH):
            newN = AsinhNorm
        elif self.smapMenu.IsChecked(ID_WATERFALL_STRETCH_SINH):
            newN = SinhNorm
        elif self.smapMenu.IsChecked(ID_WATERFALL_STRETCH_HIST):
            newN = HistEqNorm
        else:
            newN = Normalize
            
        if self.norm != newN:
            self.norm = newN
            self.draw()
            
        wx.EndBusyCursor()
        
    def onDataProduct(self, event):
        """
        Set the data product to display
        """
        
        wx.BeginBusyCursor()
        wx.Yield()
        
        for i,dmo in enumerate(self.dataMenuOptions):
            if dmo.IsChecked(ID_WATERFALL_PRODUCT_1+i):
                newIndex = i
                break
                
        if newIndex != self.index:
            self.index = newIndex
            self.draw()
            
        wx.EndBusyCursor()
        
    def onAutoDecimation(self, event):
        """
        Automatically adjust the time decimation of the plot.
        """
        
        if self.dataMenu.IsChecked(ID_WATERFALL_DECIMATION_AUTO):
            self.dadj.Enable(False)
            
            wx.BeginBusyCursor()
            wx.Yield()
        
            # Suggest a time decimation factor
            newDecFactor = max([1, int(round(self.width / 10.0 / self.parent.data.meta.dt))])
        
            if newDecFactor != self.decFactor:
                self.decFactor = newDecFactor
                self.draw()
            
            wx.EndBusyCursor()
        else:
            self.dadj.Enable(True)
            
    def onAdjustDecimation(self, event):
        """
        Manually adjust the time decimation of the plot.
        """
        
        WaterfallDecimationAdjust(self)
        
    def onBandpassOn(self, event):
        """
        Enable bandpass correction.
        """
        
        wx.BeginBusyCursor()
        wx.Yield()
        
        self.bandpass = True
        
        self.draw()
        
        wx.EndBusyCursor()
        
    def onBandpassOff(self, event):
        """
        Disable bandpass correction.
        """
        
        wx.BeginBusyCursor()
        wx.Yield()
        
        self.bandpass = False
        
        self.draw()
        
        wx.EndBusyCursor()
        
    def onShowSweep(self, event):
        """
        Toggle the pulse boundaries on/off.
        """
        
        wx.BeginBusyCursor()
        wx.Yield()
        
        self.sweep = not self.sweep
        
        self.draw()
        
        wx.EndBusyCursor()
        
    def onShowProfile(self, event):
        """
        Toggle the integrated profile on/off.
        """
        
        wx.BeginBusyCursor()
        wx.Yield()
        
        self.profile = not self.profile
        
        self.draw()
        
        wx.EndBusyCursor()
        
    def addSubplotAxes(self, fig, ax, rect, axisbg='w'):
        """
        Add a small sub-plot to an existing figure.
        
        From:
            https://stackoverflow.com/questions/17458580/embedding-small-plots-inside-subplots-in-matplotlib
        """
        
        box = ax.get_position()
        width = box.width
        height = box.height
        inax_position  = ax.transAxes.transform(rect[0:2])
        transFigure = fig.transFigure.inverted()
        infig_position = transFigure.transform(inax_position)    
        x = infig_position[0]
        y = infig_position[1]
        width *= rect[2]
        height *= rect[3]  # <= Typo was here
        subax = fig.add_axes([x,y,width,height],axisbg=axisbg)
        x_labelsize = subax.get_xticklabels()[0].get_size()
        y_labelsize = subax.get_yticklabels()[0].get_size()
        x_labelsize *= rect[2]**0.5
        y_labelsize *= rect[3]**0.5
        subax.xaxis.set_tick_params(labelsize=x_labelsize)
        subax.yaxis.set_tick_params(labelsize=y_labelsize)
        return subax
        
    def load(self):
        """
        Compute and save everything needed for the plot.
        """
        
        print("Loading PSRFITS metadata...")
        hdulist = astrofits.open(self.fitsname, mode='readonly', memmap=True)
                    
        ## File specifics
        epoch = float(hdulist[0].header['EQUINOX'])
        tStart = astro.utcjd_to_unix(hdulist[0].header['STT_IMJD'] + (hdulist[0].header['STT_SMJD'] + hdulist[0].header['STT_OFFS'])/86400.0 + astro.MJD_OFFSET)
        cFreq = hdulist[0].header['OBSFREQ']*1e6	# MHz -> Hz
        srate = hdulist[0].header['OBSBW']*1e6		# MHz -> Hz
        LFFT = hdulist[1].header['NCHAN']
        tInt = hdulist[1].header['TBIN']
        nSubs = hdulist[1].header['NSBLK']
        tSubs = nSubs*tInt
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
        freq = hdulist[1].data[0][12]*1e6
        self.freq = freq
        
        ## Pulse location
        tSweep = delay(freq, self.dm)
        self.tSweep = tSweep
        tStart = self.t
        tStop  = tStart + tSweep.max()
        subIntStart = int(tStart / tSubs) - 1
        subIntStop  = int(tStop / tSubs) + 1
        subIntStart = max([0, subIntStart])
        subIntStop  = min([subIntStop, nChunks-1])
        
        ## Spectra extraction
        print("Extracting event region...")
        samp, tRel, spec, mask = [], [], [], []
        for i in range(subIntStart, subIntStop+1):
            ### Access the correct subintegration
            subint = hdulist[1].data[i]
            
            ### Pull out various bits that we need, including:
            ###  * the start time of the subint. - tOff
            ###  * the weight mask, converted to binary - msk
            ###  * the scale and offset values - bscl and bzero
            ###  * the actual data - data
            tOff = subint[1] - subint[0]//2
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
                s = i*nSubs + j
                t = subint[1] - self.t + tInt*(j-nSubs//2)
                d = data[:,:,j]*bscl + bzero
                samp.append( s )
                tRel.append( t )
                spec.append( d )
                mask.append( msk )
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
            
        ### Come up with an appropriate smoothing window (wd) and order (od)
        ws = int(round(self.spec.shape[2]/10.0))
        ws = min([41, ws])
        if ws % 2 == 0:
            ws += 1
        od = min([9, ws-2])
        
        bpm2 = []
        for i in range(self.spec.shape[1]):
            bpm = savitzky_golay(meanSpec[i,:], ws, od, deriv=0)
            bpm = numpy.ma.array(bpm, mask=~numpy.isfinite(bpm))
            
            if bpm.mean() == 0:
                bpm += 1
            bpm2.append( bpm / bpm.mean() )
            
        ### Apply the bandpass correction
        bpm2 = numpy.array(bpm2)
        self.specBandpass = numpy.ma.array(self.spec.data*1.0, mask=self.spec.mask)
        try:
            from _helper import FastAxis0Bandpass
            FastAxis0Bandpass(self.specBandpass.data, bpm2.astype(numpy.float32))
        except ImportError:
            for i in range(self.spec.shape[1]):
                self.specBandpass.data[:,i,:] = self.spec.data[:,i,:] / bpm2[i]
                
        # Downselect to something that centers the pulse
        tPad = 0.2 + self.width
        valid = numpy.where( (self.tRel > -tPad) & (self.tRel < (tStop-tStart+tPad)) )[0]
        self.tRel = self.tRel[valid]
        self.spec = self.spec[valid,:,:]
        self.specBandpass = self.specBandpass[valid,:,:]
        
        # Run the incoherent dedispersion on the data
        print("Dedispersing data...")
        self.specD = self.spec*0
        self.specBandpassD = self.specBandpass*0
        for i in range(self.specD.shape[1]):
            self.specD[:,i,:] = incoherent(freq, self.spec[:,i,:], tInt, self.dm, boundary='fill', fill_value=numpy.nan)
            self.specBandpassD[:,i,:] = incoherent(freq, self.specBandpass[:,i,:], tInt, self.dm, boundary='fill', fill_value=numpy.nan)
        
        # Calculate the plot limits
        print("Setting default plot limits...")
        self.limits = [None,]*self.spec.shape[1]
        self.limitsBandpass = [None,]*self.spec.shape[1]
        
        try:
            from _helper import FastAxis1MinMax
            limits0 = FastAxis1MinMax(self.spec)
            limits1 = FastAxis1MinMax(self.specBandpass, chanMin=self.spec.shape[2]//10, chanMax=9*self.spec.shape[2]//10)
            if self.usedB:
                limits0 = to_dB(limits0)
                limits1 = to_dB(limits1)
            for i in range(self.spec.shape[1]):
                self.limits[i] = list(limits0[i,:])
                self.limitsBandpass[i] = list(limits1[i,:])
        except ImportError:
            toUse = range(self.spec.shape[2]//10, 9*self.spec.shape[2]//10+1)
            for i in range(self.spec.shape[1]):
                self.limits[i] = findLimits(self.spec[:,i,:], usedB=self.usedB)
            for i in range(self.spec.shape[1]):
                self.limitsBandpass[i] = findLimits(self.specBandpass[:,i,toUse], usedB=self.usedB)
                
        # Flip the axis to make the pulsar people happy
        self.spec = self.spec.T
        self.specBandpass = self.specBandpass.T
        
        # Suggest a time decimation factor
        self.decFactor = max([1, int(round(self.width / 10.0 / self.parent.data.meta.dt))])
        
        print("Ready")
        return True
        
    def render(self):
        # Clear the old figure
        self.figure.clf()
        
        self.connect()
        
    def draw(self):
        """
        Populate the figure/canvas areas with a plot.  We only need to do this
        once for this type of window.
        """
        
        # Setup
        tRel = self.tRel
        tSweep = self.tSweep
        freq = self.freq
        dataProduct = self.data_products[self.index]
        if self.bandpass:
            spec = self.specBandpass[:,self.index,:]
            specD = self.specBandpassD[:,self.index,:]
            limits = self.limitsBandpass[self.index]
        else:
            spec = self.spec[:,self.index,:]
            specD = self.specD[:,self.index,:]
            limits = self.limits[self.index]
            
        # Decimate
        nKeep = (tRel.size//self.decFactor)*self.decFactor
        tRel = tRel[:nKeep]
        tRel.shape = (nKeep//self.decFactor, self.decFactor)
        tRel = tRel.mean(axis=1)
        spec = spec[:,:nKeep]
        spec.shape = (spec.shape[0], nKeep//self.decFactor, self.decFactor)
        spec = spec.mean(axis=2)
        
        # Plot Waterfall
        self.figure.clf()
        self.ax1 = self.figure.gca()
        if self.usedB:
            m = self.ax1.imshow(to_dB(spec), interpolation='nearest', extent=(tRel[0], tRel[-1], freq[0]/1e6, freq[-1]/1e6), origin='lower', cmap=self.cmap, norm=self.norm(limits[0], limits[1]))
            try:
                cm = self.figure.colorbar(m, use_gridspec=True)
            except:
                if len(self.frame.figure1a.get_axes()) > 1:
                    self.frame.figure1a.delaxes( self.frame.figure1a.get_axes()[-1] )
                cm = self.figure.colorbar(m)
            cm.ax.set_ylabel('PSD [arb. dB]')
        else:
            m = self.ax1.imshow(spec, interpolation='nearest', extent=(tRel[0], tRel[-1], freq[0]/1e6, freq[-1]/1e6), origin='lower', cmap=self.cmap, norm=self.norm(limits[0], limits[1]))
            try:
                cm = self.figure.colorbar(m, use_gridspec=True)
            except:
                if len(self.frame.figure1a.get_axes()) > 1:
                    self.frame.figure1a.delaxes( self.frame.figure1a.get_axes()[-1] )
                cm = self.figure.colorbar(m)
            cm.ax.set_ylabel('PSD [arb. lin.]')
            
        ## Pulse boundary markers
        if self.sweep:
            self.ax1.plot(tSweep-2*self.width, freq/1e6, linestyle='--', color='r')
            self.ax1.plot(tSweep+2*self.width, freq/1e6, linestyle='--', color='r')
            
        ## Dedispersed profile
        if self.profile:
            prof = specD.sum(axis=1)
            prof = prof[:nKeep]
            prof.shape = (nKeep//self.decFactor, self.decFactor)
            prof = prof.mean(axis=1)
            ax = self.addSubplotAxes(self.figure, self.ax1, [0.7, 0.7, 0.25, 0.25])
            ax.plot(tRel, prof)
            
        self.ax1.axis('auto')
        self.ax1.set_xlim((tRel[0], tRel[-1]))
        self.ax1.set_ylim((freq[0]/1e6, freq[-1]/1e6))
        self.ax1.set_xlabel('Time - %.4f s' % self.t)
        self.ax1.set_ylabel('Frequency [MHz]')
        self.ax1.set_title(dataProduct)
        
        ## Draw and save the click (Why?)
        self.canvas.draw()
        
    def connect(self):
        """
        Connect to all the events we need to interact with the plots.
        """
        
        self.cidmotion  = self.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)
    
    def on_motion(self, event):
        """
        Deal with motion events in the stand field window.  This involves 
        setting the status bar with the current x and y coordinates as well
        as the stand number of the selected stand (if any).
        """

        if event.inaxes:
            clickX = event.xdata
            clickY = event.ydata
            
            dataX = numpy.where(numpy.abs(clickY-self.freq/1e6) == (numpy.abs(clickY-self.freq/1e6).min()))[0][0]
            dataY = numpy.where(numpy.abs(clickX-self.tRel) == (numpy.abs(clickX-self.tRel).min()))[0][0]
            
            value = to_dB(self.spec[dataX, 0, dataY])
            self.statusbar.SetStatusText("t=%.6f s, f=%.4f MHz, p=%.2f dB" % (clickX, clickY, value))
            
        else:
            self.statusbar.SetStatusText("")
    
    def disconnect(self):
        """
        Disconnect all the stored connection ids.
        """
        
        self.figure.canvas.mpl_disconnect(self.cidmotion)
        
    def onCancel(self, event):
        self.Close()
        
    def resizePlots(self, event):
        # Get the current size of the window and the navigation toolbar
        w, h = self.GetClientSize()
        wt, ht = self.toolbar.GetSize()
        
        dpi = self.figure.get_dpi()
        newW = 1.0*w/dpi
        newH = 1.0*(h-ht)/dpi
        self.figure.set_size_inches((newW, newH))
        self.figure.tight_layout()
        self.figure.canvas.draw()

    def GetToolBar(self):
        # You will need to override GetToolBar if you are using an 
        # unmanaged toolbar in your frame
        return self.toolbar


ID_WATERFALL_CONTRAST_UPR_INC = 400
ID_WATERFALL_CONTRAST_UPR_DEC = 401
ID_WATERFALL_CONTRAST_LWR_INC = 402
ID_WATERFALL_CONTRAST_LWR_DEC = 403
ID_WATERFALL_CONTRAST_OK = 404

class WaterfallContrastAdjust(wx.Frame):
    def __init__ (self, parent):	
        wx.Frame.__init__(self, parent, title='Waterfall Contrast Adjustment')
        
        self.parent = parent
        
        self.initUI()
        self.initEvents()
        self.Show()
        
        self.parent.cAdjust = self
        
    def initUI(self):
        row = 0
        panel = wx.Panel(self)
        sizer = wx.GridBagSizer(0, 0)
        
        pol = self.parent.data_products[self.parent.index]
        if self.parent.bandpass:
            typ = wx.StaticText(panel, label='%s - Bandpass' % (pol,))
        else:
            typ = wx.StaticText(panel, label='%s' % (pol,))
        sizer.Add(typ, pos=(row+0, 0), span=(1,4), flag=wx.EXPAND|wx.ALL, border=5)
        
        row += 1
        
        upr = wx.StaticText(panel, label='Upper Limit:')
        uprText = wx.TextCtrl(panel)
        uprDec = wx.Button(panel, ID_WATERFALL_CONTRAST_UPR_DEC, '-', size=(56, 28))
        uprInc = wx.Button(panel, ID_WATERFALL_CONTRAST_UPR_INC, '+', size=(56, 28))
        
        lwr = wx.StaticText(panel, label='Lower Limit:')
        lwrText = wx.TextCtrl(panel)
        lwrDec = wx.Button(panel, ID_WATERFALL_CONTRAST_LWR_DEC, '-', size=(56, 28))
        lwrInc = wx.Button(panel, ID_WATERFALL_CONTRAST_LWR_INC, '+', size=(56, 28))
        
        rng = wx.StaticText(panel, label='Range:')
        rngText = wx.TextCtrl(panel, style=wx.TE_READONLY)
        
        sizer.Add(upr,     pos=(row+0, 0), span=(1, 1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
        sizer.Add(uprText, pos=(row+0, 1), span=(1, 1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
        sizer.Add(uprDec,  pos=(row+0, 2), span=(1, 1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
        sizer.Add(uprInc,  pos=(row+0, 3), span=(1, 1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
        sizer.Add(lwr,     pos=(row+1, 0), span=(1, 1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
        sizer.Add(lwrText, pos=(row+1, 1), span=(1, 1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
        sizer.Add(lwrDec,  pos=(row+1, 2), span=(1, 1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
        sizer.Add(lwrInc,  pos=(row+1, 3), span=(1, 1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
        sizer.Add(rng,     pos=(row+2, 0), span=(1, 1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
        sizer.Add(rngText, pos=(row+2, 1), span=(1, 1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
        
        line = wx.StaticLine(panel)
        sizer.Add(line, pos=(row+3, 0), span=(1, 4), flag=wx.EXPAND|wx.BOTTOM, border=10)
        
        row += 4
        
        #
        # Buttons
        #
        
        ok = wx.Button(panel, ID_WATERFALL_CONTRAST_OK, 'Ok', size=(56, 28))
        sizer.Add(ok, pos=(row+0, 3), flag=wx.ALL, border=5)
        
        panel.SetSizer(sizer)
        sizer.Fit(self)

        self.uText = uprText
        self.lText = lwrText
        self.rText = rngText
        
        #
        # Set current values
        #
        index = self.parent.index
        if self.parent.bandpass:
            self.uText.SetValue('%.1f' % self.parent.limitsBandpass[index][1])
            self.lText.SetValue('%.1f' % self.parent.limitsBandpass[index][0])
            self.rText.SetValue('%.1f' % self.__getRange(index))
        else:
            self.uText.SetValue('%.1f' % self.parent.limits[index][1])
            self.lText.SetValue('%.1f' % self.parent.limits[index][0])
            self.rText.SetValue('%.1f' % self.__getRange(index))
        
    def initEvents(self):
        self.Bind(wx.EVT_BUTTON, self.onUpperDecrease, id=ID_WATERFALL_CONTRAST_UPR_DEC)
        self.Bind(wx.EVT_BUTTON, self.onUpperIncrease, id=ID_WATERFALL_CONTRAST_UPR_INC)
        self.Bind(wx.EVT_BUTTON, self.onLowerDecrease, id=ID_WATERFALL_CONTRAST_LWR_DEC)
        self.Bind(wx.EVT_BUTTON, self.onLowerIncrease, id=ID_WATERFALL_CONTRAST_LWR_INC)
        
        self.Bind(wx.EVT_BUTTON, self.onOk, id=ID_WATERFALL_CONTRAST_OK)
        
        self.Bind(wx.EVT_KEY_UP, self.onKeyPress)
        
    def onKeyPress(self, event):
        keycode = event.GetKeyCode()
        
        if keycode == wx.WXK_RETURN:
            index = self.parent.index
            if self.parent.bandpass:
                self.parent.limitsBandpass[index][0] = float(self.lText.GetValue())
                self.parent.limitsBandpass[index][1] = float(self.uText.GetValue())
            else:
                self.parent.limits[index][0] = float(self.lText.GetValue())
                self.parent.limits[index][1] = float(self.uText.GetValue())
            self.rText.SetValue('%.1f' % self.__getRange(index))
            self.parent.draw()
        else:
            pass
        
    def onUpperDecrease(self, event):
        index = self.parent.index
        if self.parent.bandpass:
            self.parent.limitsBandpass[index][1] -= self.__getIncrement(index)
            self.uText.SetValue('%.1f' % self.parent.limitsBandpass[index][1])
        else:
            self.parent.limits[index][1] -= self.__getIncrement(index)
            self.uText.SetValue('%.1f' % self.parent.limits[index][1])
        self.rText.SetValue('%.1f' % self.__getRange(index))
        self.parent.draw()
        
    def onUpperIncrease(self, event):
        index = self.parent.index
        if self.parent.bandpass:
            self.parent.limitsBandpass[index][1] += self.__getIncrement(index)
            self.uText.SetValue('%.1f' % self.parent.limitsBandpass[index][1])
        else:
            self.parent.limits[index][1] += self.__getIncrement(index)
            self.uText.SetValue('%.1f' % self.parent.limits[index][1])
        self.rText.SetValue('%.1f' % self.__getRange(index))
        self.parent.draw()
        
    def onLowerDecrease(self, event):
        index = self.parent.index
        if self.parent.bandpass:
            self.parent.limitsBandpass[index][0] -= self.__getIncrement(index)
            self.lText.SetValue('%.1f' % self.parent.limitsBandpass[index][0])
        else:
            self.parent.limits[index][0] -= self.__getIncrement(index)
            self.lText.SetValue('%.1f' % self.parent.limits[index][0])
        self.rText.SetValue('%.1f' % self.__getRange(index))
        self.parent.draw()
        
    def onLowerIncrease(self, event):
        index = self.parent.index
        if self.parent.bandpass:
            self.parent.limitsBandpass[index][0] += self.__getIncrement(index)
            self.lText.SetValue('%.1f' % self.parent.limitsBandpass[index][0])
        else:
            self.parent.limits[index][0] += self.__getIncrement(index)
            self.lText.SetValue('%.1f' % self.parent.limits[index][0])
        self.rText.SetValue('%.1f' % self.__getRange(index))
        self.parent.draw()
        
    def onOk(self, event):
        needToRedraw = False
        
        index = self.parent.index
        if self.parent.bandpass:
            if float(self.lText.GetValue()) != self.parent.limitsBandpass[index][0]:
                self.parent.limitsBandpass[index][0] = float(self.lText.GetValue())
                needToRedraw = True
            if float(self.uText.GetValue()) != self.parent.limitsBandpass[index][1]:
                self.parent.limitsBandpass[index][1] = float(self.uText.GetValue())
                needToRedraw = True
                
        else:
            if float(self.lText.GetValue()) != self.parent.limits[index][0]:
                self.parent.limits[index][0] = float(self.lText.GetValue())
                needToRedraw = True
            if float(self.uText.GetValue()) != self.parent.limits[index][1]:
                self.parent.limits[index][1] = float(self.uText.GetValue())
                needToRedraw = True
                
        if needToRedraw:
            self.parent.draw()
            
        self.parent.cAdjust = None
        self.Close()
        
    def __getRange(self, index):
        if self.parent.bandpass:
            return (self.parent.limitsBandpass[index][1] - self.parent.limitsBandpass[index][0])
        else:
            return (self.parent.limits[index][1] - self.parent.limits[index][0])
        
    def __getIncrement(self, index):
        return 0.1*self.__getRange(index)


ID_WATERFALL_DECIMATION_INC = 500
ID_WATERFALL_DECIMATION_DEC = 501
ID_WATERFALL_DECIMATION_OK = 502

class WaterfallDecimationAdjust(wx.Frame):
    def __init__ (self, parent):	
        wx.Frame.__init__(self, parent, title='Waterfall Decimation Adjustment')
        
        self.parent = parent
        
        self.initUI()
        self.initEvents()
        self.Show()
        
    def initUI(self):
        row = 0
        panel = wx.Panel(self)
        sizer = wx.GridBagSizer(0, 0)
        
        upr = wx.StaticText(panel, label='Time Decimation:')
        uprText = wx.TextCtrl(panel)
        uprDec = wx.Button(panel, ID_WATERFALL_DECIMATION_DEC, '-', size=(56, 28))
        uprInc = wx.Button(panel, ID_WATERFALL_DECIMATION_INC, '+', size=(56, 28))
        
        rng = wx.StaticText(panel, label='Bins/Pulse:')
        rngText = wx.TextCtrl(panel, style=wx.TE_READONLY)
        
        sizer.Add(upr,     pos=(row+0, 0), span=(1, 1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
        sizer.Add(uprText, pos=(row+0, 1), span=(1, 1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
        sizer.Add(uprDec,  pos=(row+0, 2), span=(1, 1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
        sizer.Add(uprInc,  pos=(row+0, 3), span=(1, 1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
        sizer.Add(rng,     pos=(row+1, 0), span=(1, 1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
        sizer.Add(rngText, pos=(row+1, 1), span=(1, 1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
        
        line = wx.StaticLine(panel)
        sizer.Add(line, pos=(row+2, 0), span=(1, 4), flag=wx.EXPAND|wx.BOTTOM, border=10)
        
        row += 3
        
        #
        # Buttons
        #
        
        ok = wx.Button(panel, ID_WATERFALL_DECIMATION_OK, 'Ok', size=(56, 28))
        sizer.Add(ok, pos=(row+0, 3), flag=wx.ALL, border=5)
        
        panel.SetSizer(sizer)
        sizer.Fit(self)

        self.uText = uprText
        self.rText = rngText
        
        #
        # Set current values
        #
        self.uText.SetValue('%i' % self.parent.decFactor)
        self.rText.SetValue('%.2f' % self.__getBinsPerPulse())
        
    def initEvents(self):
        self.Bind(wx.EVT_BUTTON, self.onUpperDecrease, id=ID_WATERFALL_DECIMATION_DEC)
        self.Bind(wx.EVT_BUTTON, self.onUpperIncrease, id=ID_WATERFALL_DECIMATION_INC)
        
        self.Bind(wx.EVT_BUTTON, self.onOk, id=ID_WATERFALL_DECIMATION_OK)
        
        self.Bind(wx.EVT_KEY_UP, self.onKeyPress)
        
    def onKeyPress(self, event):
        keycode = event.GetKeyCode()
        
        if keycode == wx.WXK_RETURN:
            self.parent.decFactor = int(self.uText.GetValue())
            self.rText.SetValue('%.2f' % self.__getBinsPerPulse())
            self.parent.draw()
        else:
            pass
        
    def onUpperDecrease(self, event):
        self.parent.decFactor -= 1
        if self.parent.decFactor < 1:
            self.parent.decFactor = 1
        self.uText.SetValue('%i' % self.parent.decFactor)
        self.rText.SetValue('%.2f' % self.__getBinsPerPulse())
        
        self.parent.draw()
        
    def onUpperIncrease(self, event):
        self.parent.decFactor += 1
        self.uText.SetValue('%i' % self.parent.decFactor)
        self.rText.SetValue('%.2f' % self.__getBinsPerPulse())
        
        self.parent.draw()
        
    def onOk(self, event):
        needToRedraw = False
        
        if self.parent.decFactor != int(self.uText.GetValue()):
            self.parent.decFactor = int(self.uText.GetValue())
            needToRedraw = True
            
        if needToRedraw:
            self.parent.draw()
            
        self.Close()
        
    def __getBinsPerPulse(self):
        return self.parent.width // (self.parent.decFactor * self.parent.parent.data.meta.dt)


class HtmlWindow(wx.html.HtmlWindow): 
    def __init__(self, parent): 
        wx.html.HtmlWindow.__init__(self, parent, style=wx.NO_FULL_REPAINT_ON_RESIZE|wx.SUNKEN_BORDER) 
        
        if "gtk2" in wx.PlatformInfo: 
            self.SetStandardFonts()
            
    def OnLinkClicked(self, link): 
        a = link.GetHref()
        if a.startswith('#'): 
            wx.html.HtmlWindow.OnLinkClicked(self, link) 
        else: 
            wx.LaunchDefaultBrowser(link.GetHref())


class HelpWindow(wx.Frame):
    def __init__(self, parent):
        wx.Frame.__init__(self, parent, -1, 'plotHDF Handbook', size=(570, 400))
        
        self.initUI()
        self.Show()
        
    def initUI(self):
        panel = wx.Panel(self, -1, style=wx.BORDER_SUNKEN)
        vbox = wx.BoxSizer(wx.VERTICAL)
        
        help = HtmlWindow(panel)
        #help.LoadPage(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'docs/help.html'))
        help.SetPage("""<html>

<body>
<a name="top"><h4>Table of Contents</h4></a>
<ul>
<li><a href="#intro">Introduction</a></li>
<li><a href="#layout">Window Layout</a></li>
<li><a href="#usage">Usage</a></li>
<li><a href="#mouse">Mouse Interaction</a></li>
<li><a href="#keyboard">Keyboard Interaction</a></li>
</ul>

<p>
<a name="intro">
<h6>Introduction</h6>
plotSinglePulse is a graphical interface for working with single pulse search results from 
PRESTO.  This script allows the user to plot the candidates in an interactive fashion similar
to the postscript plots generated by single_pulse_search.py as well as to plot the associated
waterfall data if a PSRFITS file is provided.
<br /><a href="#top">Top</a>
</a>
</p>

<p>
<a name="layout">
<h6>Window Layout</h6>
The plotSinglePulse window is broken into two vertical sections.  The top section shows, from 
left to right, the histogram of signal-to-noise ratio (S/N) for all pulses displayed, the 
histogram of dispersion measure (DM) for all pulses displayed, and a scatter plot of DM vs. S/N.
The lower section shows a scatter plot of time vs. DM with user-configurable colors and symbols.
<br /><br />
<b>Note:</b>  Since there may be may candidates the upper left and lower plot windows use internal
decimation to display the data.  This default behavior limits the number of points plotted in 
each window and the decimation can be adjusted via the <a href="#usage">Data menu</a>.
<br /><a href="#top">Top</a>
</a>
</p>

<p>
<a name="usage">
<h6>Usage</h6>
After the single pulse candidates are loaded file has been loaded there are a variety of menu, 
<a href="#mouse">mouse</a>, and <a href="#keyboard">keyboard commands</a> that can be used to 
interact with the data.  The menus are:
<ul>
    <li>Color - Adjust the parameter used for color coding, color stretch, map, and transfer function used in the dynamic waterfall</li>
    <li>Data - Change the symbol coding, size, and pulse selection thresholds</li>
    <li>Display - Toggle display decimation on/off and change the decimation factors</li>
    <li>Help - Show this help message.</li>
</ul>
<br /><a href="#top">Top</a>
</a>
</p>

<p>
<a name="mouse">
<h6>Mouse Interaction</h6>
The mouse can be used in any of the plotting panels.  In the upper left and center sections:
<ul>
    <li>Left Click - Nothing</li>
    <li>Middle Click - Select pulses in the current bin for highlighting in the lower panel</li>
    <li>Right Clock - Clear the last selection make with the middle click</li>
</ul>
<br /><br />
In the lower section:
<ul>
    <li>Left Click - Move the cross hairs around the screen</li>
    <li>Middle Click - Unmask the closest pulse candidate to the cross hairs</li>
    <li>Right Click - Mask the closest pulse candidate to the cross hairs</li>
</ul>
<br /><br />
The lower section can also be zoomed by selecting the zoom tool from the toolbar in the lower left 
of the window.  <b>Note:</b> While the zoom tool is active the standard mouse functions in this window
are disabled.
<br /><a href="#top">Top</a>
</a>
</p>

<p>
<a name="keyboard">
<h6>Keyboard Interaction</h6>
The keyboard can also be used to interact with the lower plotting panel.  The keyboard commands 
are:
<ul>
    <li>p - print(the information about the underlying pulse</li>)
    <li> e - export the current pulses to a file</li>
    <li>s - display a DM time slice</li>
    <li>w - display the PSRFITS waterfall for a pulse</li>
    <li>u - unmask pulses in a region of time</li>
    <li>m - mask pulses in a region of time</li>
    <li>y - unmask pulses in a region of time/DM</li>
    <li>n - mask pulses in a region of time/DM</li>
    <li>t - unmask pulses in a region of DM</li>
    <li>b - mask pulses in a region of DM</li>
</ul>
<br /><br />
<b>Note:</b> In order to interact with the lower section you will need to click in the axes with the left mouse button.
<br /><a href="#top">Top</a>
</a>
</p>

</body>

</html>""")
        vbox.Add(help, 1, wx.EXPAND)
        
        self.CreateStatusBar()
        
        panel.SetSizerAndFit(vbox)


def main(args):
    # Turn off all NumPy warnings to keep stdout clean
    errs = numpy.geterr()
    numpy.seterr(all='ignore', invalid='ignore', divide='ignore')
    
    # Parse the command line options
    args.time_range  = [float(v) for v in args.time_range.split(',') ]
    args.dm_range    = [float(v) for v in args.dm_range.split(',')   ]
    args.width_range = [float(v) for v in args.width_range.split(',')]
    
    # Check for the _helper module
    try:
        import _helper
    except ImportError:
        print("WARNING: _helper.so not found, consider building it with 'make'")
        
    # Go!
    app = wx.App(0)
    frame = MainWindow(None, -1)
    frame.data = SinglePulse_GUI(frame)
    frame.render()
    if args.filename is not None:
        ## If there is a filename on the command line, load it
        frame.filenames = args.filename
        frame.data.loadData(args.filename, threshold=args.threshold, timeRange=args.time_range, dmRange=args.dm_range, widthRange=args.width_range, fitsname=args.fitsname)
        frame.data.render()
        frame.data.draw()
        
    else:
        pass
        
    app.MainLoop()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='read in a collection of .singlepulse files and plot them interactively', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('filename', type=str, nargs='+', 
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
    
