"""
Unit tests for the a small DRX file.
"""

# Python3 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info > (3,):
    xrange = range
    
import unittest
import os
import re
import glob
import numpy
import subprocess


_URL = 'https://lda10g.alliance.unm.edu/tutorial/UnknownPulsar/056227_000024985_DRX.dat'
_SIZE_MB = 250
_FILENAME = 'data/drx.dat'


currentDir = os.path.abspath(os.getcwd())
if os.path.exists(os.path.join(currentDir, 'test_drx.py')):
    MODULE_BUILD = currentDir
else:
    MODULE_BUILD = None
    
run_scripts_tests = False
if MODULE_BUILD is not None:
    run_scripts_tests = True


class pulsar_tests(unittest.TestCase):
    def setUp(self):
        """Make sure we have the comparison files in place."""
        
        # Raw data
        if not os.path.exists(_FILENAME):
            subprocess.check_call(['curl', _URL, 
                                   '--range', '0-%i' % (int(_SIZE_MB)*1024*1024), 
                                   '-o', 'data/drx.dat', '--create-dirs'])
            
    def tearDown(self):
        for filename in glob.glob('*.fits'):
            try:
                os.unlink(filename)
            except OSError:
                pass
                
        try:
            os.unlink('script.log')
        except OSError:
            pass

def _test_generator(script):
    """
    Function to build a test method for each script that is provided.  
    Returns a function that is suitable as a method inside a unittest.TestCase
    class
    """
    
    def test(self):
        with open('script.log', 'w') as logfile:
            try:
                cmd = ['python',]
                cmd.extend(script.split())
                cmd.append(_FILENAME)
                status = subprocess.check_call(cmd, stdout=logfile)
            except subprocess.CalledProcessError:
                status = 1
                
        if status == 1:
            with open('logfile', 'r') as logfile:
                print(logfile.read())
        self.assertEqual(status, 0)
        
    return test


def _name_to_name(filename):
    filename = filename.split()[0]
    filename = os.path.splitext(filename)[0]
    parts = filename.split(os.path.sep)
    start = parts.index('..')
    parts = parts[start+1:]
    return '_'.join(parts)


if run_scripts_tests:
    _SCRIPTS = ['../drx2dri.py', 
                '../writePsrfits2.py --ra=00:00:00 --dec=00:00:00',
                '../writePsrfits2D.py --ra=00:00:00 --dec=00:00:00 4.0']
    _SCRIPTS.sort()
    for script in _SCRIPTS:
        test = _test_generator(script)
        name = 'test_%s' % _name_to_name(script)
        doc = """Simple execution of the '%s' script.""" % os.path.basename(script)
        setattr(test, '__doc__', doc)
        setattr(pulsar_tests, name, test)


class pulsar_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the DRX pulsar script
    tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(pulsar_tests))


if __name__ == '__main__':
    unittest.main()
    