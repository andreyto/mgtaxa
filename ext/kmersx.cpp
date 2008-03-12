#include "mgtaxa/kmers.hpp"

#include <boost/python.hpp>

namespace bp = boost::python;

namespace MGT {

class KmerCounterPy : public KmerCounter {

	KmerCounterPy(int kmerLen) :
	KmerCounter(kmerLen)
	{}
	
	void process(bp::array seq) {
	}

}; // class KmerCounterPy

} // namespace MGT

using namespace boost::python;

BOOST_PYTHON_MODULE(kmersx)
{
    class_<MGT::KmerCounterPy>("KmerCounter", init<int> >())
        .def("process", &MGT::KmerCounterPy::process)
        .def("getCounts", &MGT::KmerCounterPy::getCounts);
}
