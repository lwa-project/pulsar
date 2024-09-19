"""
Unit tests for the various pulsar scripts.
"""

import unittest
import glob
import sys
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
    from pylint.reporters.json_reporter import JSONReporter
    if MODULE_BUILD is not None:
        run_scripts_tests = True
except ImportError:
    pass


__version__  = "0.1"
__author__   = "Jayce Dowell"


_LINT_RE = re.compile('(?P<module>.*?):(?P<line>\d+): (error )?[\[\(](?P<type>.*?)[\]\)] (?P<info>.*)')


_PYLINT_IGNORES = [ ]


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
    """A unittest.TestCase collection of unit tests for the LSL scripts."""
    
    def test_scripts(self):
        """Static analysis of the LSL scripts."""
        
        _SCRIPTS = glob.glob(os.path.join(MODULE_BUILD, '..', 'scripts', '*.py'))
        _SCRIPTS.sort()
        for script in _SCRIPTS:
            name = script.rsplit('scripts'+os.path.sep)[-1]
            with self.subTest(script=name):
                pylint_output = StringIO()
                reporter = JSONReporter(pylint_output)
                pylint_args = [script, "-E", "--extension-pkg-whitelist=numpy,ephem,lsl", "--init-hook='import sys; sys.path=[%s]; sys.path.insert(0, \"%s\")'" % (",".join(['"%s"' % p for p in sys.path]), os.path.dirname(MODULE_BUILD))]
                Run(pylint_args, reporter=reporter, exit=False)
                results = json.loads(pylint_output.getvalue())
                
                for i,entry in enumerate(results):
                    with self.subTest(error_number=i+1):
                        false_positive = False
                        for isym,imesg in _PYLINT_IGNORES:
                            if entry['symbol'] == isym and entry['message'].startswith(imesg):
                                false_positive = True
                                break
                        if false_positive:
                            continue
                            
                        self.assertTrue(False, f"{entry['path']}:{entry['line']} - {entry['symbol']} - {entry['message']}")


class scripts_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the pulsar script
    tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(scripts_tests))


if __name__ == '__main__':
    unittest.main()
    
