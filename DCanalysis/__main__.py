import argparse
import unittest

parser = argparse.ArgumentParser()
parser.add_argument('-test', '--test', action='store_true', default=False,
                    help='run unittests')
args = parser.parse_args()


if args.test:
    testsuite = unittest.TestLoader().discover('.')
    unittest.TextTestRunner(verbosity=1).run(testsuite)
