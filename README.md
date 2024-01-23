[![GHA](https://github.com/lwa-project/pulsar/actions/workflows/main.yml/badge.svg)](https://github.com/lwa-project/pulsar/actions/workflows/main.yml)

[![Paper](https://img.shields.io/badge/arXiv-1410.7422-blue.svg)](https://arxiv.org/pdf/1410.7422.pdf)

LWA1 Pulsar Scripts
===================
The LSL Pulsar extension contains a variety of tools for converting LWA1 
data to the PSRFITS format.  These tools build off the LSL 0.6.x/1.0.x 
frameworks as well as the psrfits_utils Python module.  psrfits_utils can
be found at:  https://github.com/lwa-project/psrfits_utils

No installation (e.g., python setup.py install) is required to use the 
software.  Simply run make to build the `_psr.so` and `_helper.so` modules and 
then use the scripts in the Pulsar directory.

writePsrfits2.py
----------------
Given a DRX file, create a PSRFITS file of the data.

writePsrfits2D.py
-----------------
Given a DRX file, coherently dedisperse the data at the specified DM and 
create a PSRFITS file of the data.

writePsrfits2FromDRSpec.py
--------------------------
Given a DR spectrometer file in either XX/YY, IV, or IQUV mode, create a 
PSRFITS file of the data.

writePsrfits2FromHDF5.py
------------------------
Given an HDF5 created by hdfWaterfall.py or drspec2hdf.py, create a PSRFITS 
file of the data.

writePsrfits2Multi.py
---------------------
Experimental script to that takes in a collection of DRX files observed at 
the same time and processes them such that they are aligned in time.  This
yields PSRFITS files that can be combined with the 'combine_lwa' script
across multiple beams.

writePsrfits2DMulti.py
----------------------
Experimental script to that takes in a collection of DRX files observed at 
the same time and processes them such that they are both aligned in time and
coherently dedispersed.  This yields PSRFITS files that can be combined with
the 'combine_lwa' script across multiple beams.

writeHDF5FromPsrfits.py
-----------------------
Given a PSRFITS file created by the writePsrfit2 family of converters, build
an HDF5 file that is compatible with the tools in Commissioning/DRX/HDF5.

updatePsrfitsMask.py
--------------------
Use spectral kurtosis to update the weight mask in a PSRFITS file to flag RFI.
The script also takes in a list of frequencies/frequency ranges that can be 
used to update the weight mask.

plotSinglePulse.py
------------------
Graphical interface for working with .singlepulse search results from PRESTO.

dedispersion.cpp/fft.cpp/kurtosis.cpp/psr.cpp/quantize.cpp/reduce.cpp/utils.cpp
-------------------------------------------------------------------------------
Compiled Python/NumPy extension used by the scripts to speed up processing.

helper.cpp
----------
Compiled Python/NumPy extension used by the script to speed up display.

Makefile
--------
Makefile for the psr.c module.

setup.py
--------
File used by the Makefile to help build the Python extensions.

data.py
-------
External file from Commissioning/DRX/HDF5 used by writeHDF5FromPsrfits.py
for building the HDF5 files.

drx2drxi.py
-----------
Script to take a standard DRX file and convert it to two polarization interlaced
DRX files, one for each tuning.  The interlaced DRX format is compatible with
the dspsr suite.
