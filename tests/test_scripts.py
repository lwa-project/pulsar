"""
Unit tests for the various pulsar scripts.
"""

import unittest
import glob
import sys
import imp
import re
import os

currentDir = os.path.abspath(os.getcwd())
if os.path.exists(os.path.join(currentDir, 'test_scripts.py')):
    MODULE_BUILD = os.path.join(currentDir, '..')
else:
    MODULE_BUILD = None
    
run_scripts_tests = False
try:
    from io import StringIO
    from pylint.lint import Run
    from pylint.reporters.text import TextReporter
    if MODULE_BUILD is not None:
        run_scripts_tests = True
except ImportError:
    pass


__version__  = "0.1"
__author__   = "Jayce Dowell"


_LINT_RE = re.compile('(?P<module>.*?)\:(?P<line>\d+)\: (error )?[\[\(](?P<type>.*?)[\]\)] (?P<info>.*)')


_SAFE_TO_IGNORE = ["Possible",
                   "Module 'numpy",
                   "Module 'ephem",
                   "Module 'wx",
                   "Unable to import 'wx",
                   "Unable to import 'presto",
                   "No name 'erf' in module 'scipy.special'",
                   "Instance of 'HDUList' has no 'header' member",
                   "Instance of 'HDUList' has no 'data' member",
                   "Instance of 'Group' has no 'shape' member",
                   "Undefined variable 'BindToCore",
                   "Undefined variable 'BindOpenMPToCores",
                   "Undefined variable 'PulsarEngineRaw",
                   "Undefined variable 'PulsarEngineRawWindow",
                   "Undefined variable 'PhaseRotator",
                   "Undefined variable 'ComputeSKMask",
                   "Undefined variable 'ComputePseudoSKMask",
                   "Undefined variable 'MultiChannelCD",
                   "Undefined variable 'CombineToIntensity",
                   "Undefined variable 'CombineToLinear",
                   "Undefined variable 'CombineToCircular",
                   "Undefined variable 'CombineToStokes",
                   "Undefined variable 'OptimizeDataLevels8Bit",
                   "Undefined variable 'OptimizeDataLevels4Bit",
                   "Undefined variable 'useWisdom",
                   "Class 'int' has no 'from_bytes' member",
                   "Module 'astropy.units' has no 'hourangle' member",
                   "Module 'astropy.units' has no 'degree' member",
                   "Instance of 'Empty' has no 'decode' member"]


def _get_context(filename, line, before=0, after=0):
    to_save = range(line-1-before, line-1+after+1)
    context = []
    with open(filename, 'r') as fh:
        i = 0
        for line in fh:
            if i in to_save:
                context.append(line)
            i += 1
    return context


@unittest.skipUnless(run_scripts_tests, "requires the 'pylint' module")
class scripts_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the pulsar scripts."""
    
    pass     


def _test_generator(script):
    """
    Function to build a test method for each script that is provided.  
    Returns a function that is suitable as a method inside a unittest.TestCase
    class
    """
    
    def test(self):
        pylint_output = StringIO()
        reporter = TextReporter(pylint_output)
        pylint_args = [script, "-E", "--extension-pkg-whitelist=numpy,ephem,lsl", "--init-hook='import sys; sys.path=[%s]; sys.path.insert(0, \"%s\")'" % (",".join(['"%s"' % p for p in sys.path]), os.path.dirname(MODULE_BUILD))]
        Run(pylint_args, reporter=reporter)
        out = pylint_output.getvalue()
        out_lines = out.split('\n')
        
        for line in out_lines:
            ignore = False
            for phrase in _SAFE_TO_IGNORE:
                if line.find(phrase) != -1:
                    ignore = True
                    break
            if ignore:
                continue
                
            mtch = _LINT_RE.match(line)
            if mtch is not None:
                line_no, type, info = mtch.group('line'), mtch.group('type'), mtch.group('info')
                ignore = False
                if line.find('before assignment') != -1:
                    context = _get_context(script, int(line_no), before=20, after=20)
                    loc = context[20]
                    level = len(loc) - len(loc.lstrip()) - 4
                    if loc.strip().rstrip() == 'rawSpectraDedispersed)':
                        ## This is fragile on purpose
                        loc = context[20-4]
                        level = len(loc) - len(loc.lstrip()) - 4
                    found_try = None
                    for i in range(1, 20):
                        if context[20-i][level:level+3] == 'try':
                            found_try = i
                            break
                    found_except = None
                    for i in range(0, 20):
                        if context[20+i][level:level+6] == 'except':
                            found_except = i
                            break
                    if found_try is not None and found_except is not None:
                        if context[20+found_except].find('NameError') != -1:
                            ignore = True
                if ignore:
                    continue
                    
                self.assertEqual(type, None, "%s:%s - %s" % (os.path.basename(script), line_no, info))
    return test


if run_scripts_tests:
    _SCRIPTS = glob.glob(os.path.join(MODULE_BUILD, '*.py'))
    for depth in range(1, 3):
        path = [MODULE_BUILD, '..']
        path.extend(['*',]*depth)
        path.append('*.py')
        _SCRIPTS.extend(glob.glob(os.path.join(*path)))
    _SCRIPTS = list(filter(lambda x: x.find('test_scripts.py') == -1, _SCRIPTS))
    _SCRIPTS = list(filter(lambda x: x.find('psrfits_utils') == -1, _SCRIPTS))
    _SCRIPTS.sort()
    for script in _SCRIPTS:
        test = _test_generator(script)
        name = 'test_%s' % os.path.splitext(os.path.basename(script))[0]
        doc = """Static analysis of the '%s' script.""" % os.path.basename(script)
        setattr(test, '__doc__', doc)
        setattr(scripts_tests, name, test)


class scripts_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the pulsar script
    tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(scripts_tests))


if __name__ == '__main__':
    unittest.main()
    
