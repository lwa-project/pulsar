"""
Unit tests for a small HDF5 files generated from a DR spectrometer file.
"""

import unittest
import os
import re
import sys
import glob
import numpy
import subprocess
from urllib.request import Request, urlopen


_URL = 'https://lda10g.alliance.unm.edu/tutorial/B0329+54/056770_000044687'
_SIZE_MB = 250
_FILENAME = 'data/drspec-waterfall.hdf5'


currentDir = os.path.abspath(os.getcwd())
if os.path.exists(os.path.join(currentDir, 'test_pulsar_hdf5.py')):
    MODULE_BUILD = currentDir
else:
    MODULE_BUILD = None
    
run_scripts_tests = False
if MODULE_BUILD is not None:
    run_scripts_tests = True


class pulsar_hdf5_tests(unittest.TestCase):
    def setUp(self):
        """Make sure we have the comparison files in place."""
        
        # Raw data
        if not os.path.exists(_FILENAME):
            os.makedirs(os.path.dirname(_FILENAME), exist_ok=True)
            
            req = Request(_URL)
            req.add_header("Range", f"bytes=0-{_SIZE_MB*1024*1024}")
            with open(os.path.join(os.path.dirname(_FILENAME), 'drspec.dat'), 'wb') as fh:
                with urlopen(req) as uh:
                    while True:
                        data = uh.read(16*1024*1024)
                        if data:
                            fh.write(data)
                        else:
                            break
                            
            req = Request('https://raw.githubusercontent.com/lwa-project/commissioning/main/DRX/HDF5/drspec2hdf.py')
            with open(os.path.join(os.path.dirname(_FILENAME), 'drspec2hdf.py'), 'wb') as fh:
                with urlopen(req) as uh:
                    while True:
                        data = uh.read(16*1024*1024)
                        if data:
                            fh.write(data)
                        else:
                            break
            subprocess.check_call(['cp', '../data.py', os.path.join(os.path.dirname(_FILENAME), 'data.py')])
            subprocess.check_call([sys.executable, 'drspec2hdf.py', 'drspec.dat'],
                                   cwd=os.path.dirname(_FILENAME))
            
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
                cmd = [sys.executable,]
                cmd.extend(script.split())
                cmd.append(_FILENAME)
                status = subprocess.check_call(cmd, stdout=logfile)
            except subprocess.CalledProcessError:
                status = 1
                
        if status == 1:
            with open('script.log', 'r') as logfile:
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
    _SCRIPTS = ['../writePsrfits2FromHDF5.py --source=B0329+54 --ra=3:32:59.368 --dec=54:34:43.57',]
    _SCRIPTS.sort()
    for script in _SCRIPTS:
        test = _test_generator(script)
        name = 'test_%s' % _name_to_name(script)
        doc = """Simple execution of the '%s' script.""" % os.path.basename(script.split()[0])
        setattr(test, '__doc__', doc)
        setattr(pulsar_hdf5_tests, name, test)


class pulsar_hdf5_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the HDF5
    pulsar script tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(pulsar_hdf5_tests))


if __name__ == '__main__':
    unittest.main()
