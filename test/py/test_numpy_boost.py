# Copyright David Abrahams 2004. Distributed under the Boost
# Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

false = 0;
true = 1;

import doctest, boost_numpy_tests
def _count_failures(test_modules = (boost_numpy_tests,)):
    failures = 0
    for m in test_modules:
        failures += doctest.testmod(m)[0]
    return failures

def _run(args = None):
    import sys #, boost_numpy_tests

    if args is not None:
        sys.argv = args

    import numpy
    m = numpy

    # test the info routine outside the doctest. See numpy.cpp for an
    # explanation
    import test_numpy_boostx
    test_numpy_boostx.set_module_and_type('numpy', 'ndarray')
    #test_numpy_boostx.info(m.array((1,2,3)))
    a = test_numpy_boostx.new_array()
    #print a, m.info(a)
    b = m.ones((3,3),float)
    print b, m.info(b)
    return
    from boost_printer import printer
    from test_numpy_boostx import *
    x = new_array()
    y = x.copy()
    p = printer()
    check = p.check
    exercise_numarray(x, p)	

    failures = 0

    #
    # Run tests 4 different ways if both modules are installed, just
    # to show that set_module_and_type() is working properly
    #
    
    # run all the tests with default module search
    print 'testing default extension module:', \
          test_numpy_boostx.get_module_name() or '[numeric support not installed]'

    failures += _count_failures()
        
    print 'testing Numpy module explicitly'
    test_numpy_boostx.set_module_and_type('numpy', 'ndarray')
        
    failures += _count_failures()
            

    # see that we can go back to the default
    test_numpy_boostx.set_module_and_type('', '')
    print 'testing default module again:', \
          test_numpy_boostx.get_module_name() or '[numeric support not installed]'
    
    failures += _count_failures()
    
    return failures
    
if __name__ == '__main__':
    print "running..."
    import sys
    status = _run()
    if (status == 0): print "Done."
    sys.exit(status)
