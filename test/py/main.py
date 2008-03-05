import unittest

# list all test modules (without .py extension)
tests = [
    'test_sample_boost',
	]

# import all test modules
for t in tests:
    exec 'import ' + t


def main():
    # load all tests suites
    suites = [ eval(t + '.suite()') for t in tests ]

    # and make global test suite
    ts = unittest.TestSuite( suites )

    # finally run it
    unittest.TextTestRunner().run( ts )

if __name__ == '__main__':
    main()