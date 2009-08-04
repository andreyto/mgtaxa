//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the MGTAXA package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
char const* greet()
{
   return "hello, world";
}

#include <boost/python.hpp>

//If python.h has to be included, it has to be done through this:
//#include "boost/python/detail/wrap_python.hpp"
//and before any before any other system headers
//(from http://boost.org/libs/python/doc/building.html#configuring-boost-build)

using namespace boost::python;

BOOST_PYTHON_MODULE(test_sample_boostx)
{
    def("greet", greet);
}
