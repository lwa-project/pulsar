"""
Unit tests for the various pulsar scripts.
"""

import unittest
import glob
import sys
import os
import json

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
        
        # Pre-seed data.py
        os.system("%s ../data.py" % sys.executable)
except ImportError:
    pass


__version__  = "0.2"
__author__   = "Jayce Dowell"


_PYLINT_IGNORES = [('no-member',    "Instance of 'HDUList'"),
                   ('no-member',    "Instance of 'Group' has no"),
                   ('no-member',    "Module 'wx' has no"),
                   ('no-member',    "Module 'wx.html' has no"),
                   ('import-error', "Unable to import 'presto")]


@unittest.skipUnless(run_scripts_tests, "requires the 'pylint' module")
class scripts_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the LSL scripts."""
    
    def test_scripts(self):
        """Static analysis of the pulsar scripts."""
        
        _SCRIPTS = glob.glob(os.path.join(MODULE_BUILD, '*.py'))
        _SCRIPTS.sort()
        for script in _SCRIPTS:
            name = script.rsplit('scripts'+os.path.sep)[-1]
            with self.subTest(script=name):
                pylint_output = StringIO()
                reporter = JSONReporter(pylint_output)
                pylint_args = [script, "-E", "--extension-pkg-allow-list=_psr", "--init-hook='import sys; sys.path=[%s]; sys.path.insert(0, \"%s\")'" % (",".join(['"%s"' % p for p in sys.path]), os.path.dirname(MODULE_BUILD))]
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
    
