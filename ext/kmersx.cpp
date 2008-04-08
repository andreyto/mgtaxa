#include "mgtaxa/kmers.hpp"

#define PY_ARRAY_UNIQUE_SYMBOL PyArrayHandle

#include <boost/python.hpp>

#include "mgtaxa/py_num_util.hpp"

namespace bp = boost::python;
namespace bpn = bp::numeric;

namespace MGT {

/** Obtain data pointer cast to a specific datatype, size and stride for a 1d Numpy array in one call.
*/

template<typename CType>
char *
getCData1D(bpn::array a,NPY_TYPES npyType, npy_intp& size, npy_intp& stride) {
    num_util::check_rank(a,1);
    //num_util::check_type(a,npyType);
    if( sizeof(CType) != PyArray_ITEMSIZE(a.ptr()) ){
        PyErr_SetString(PyExc_ValueError, "C type size does not match NumPy array item size");
        bp::throw_error_already_set();
    }
    size = num_util::size(a);
    stride = PyArray_STRIDE(a.ptr(),0);
    //ADT_LOG << ADT_OUTVAR(size) << ADT_OUTVAR(stride) << '\n';
    return reinterpret_cast<char*>(num_util::data(a));
}

class KmerCounterPy : public KmerCounter {

public:
	
    KmerCounterPy(int kmerLen) :
	KmerCounter(kmerLen),
    m_seqLen(0)
	{
        //ADT_LOG << ADT_OUTVAR(kmerLen) << '\n';
    }
	
	void process(bpn::array seq) {
        npy_intp nSeq = 0, stride = 0;
        const char *pSeq = getCData1D<char>(seq,NPY_CHAR,nSeq,stride);
        for(npy_intp i = 0; i < nSeq; i++) {
            doCNuc(pSeq[i*stride]);
        }
        m_seqLen += nSeq;
	}

    bp::tuple counts(bpn::array counts, bpn::array indices) {
        //ADT_LOG << ADT_OUTVAR(m_seqLen) << '\n';
        npy_intp totalCount = 0;
        npy_intp sizeCounts = 0, sizeIndices = 0, strideCounts = 0, strideIndices = 0;
        char * pCounts = getCData1D<npy_int32>(counts,NPY_INT32,sizeCounts,strideCounts);
        char * pIndices = getCData1D<npy_int64>(indices,NPY_INT64,sizeIndices,strideIndices);
        if( m_seqLen >  sizeCounts || m_seqLen > sizeIndices){
            PyErr_SetString(PyExc_ValueError, "The length of output array should be no less than the length of processed sequence");
            bp::throw_error_already_set();
        }
        npy_intp outSize = numKmers();
        //ADT_LOG << ADT_OUTVAR(numKmers()) << '\n';
        if( outSize >  sizeCounts ){
            PyErr_SetString(PyExc_ValueError, "Output size execeeds array bounds");
            bp::throw_error_already_set();
        }
        startKmer();
        for(npy_intp i = 0; i < outSize; i++,nextKmer()) {
            int cnt = getKmerCount();
            int id = getKmerId();
            (*reinterpret_cast<npy_int64*>(pIndices+i*strideIndices)) = id;
            (*reinterpret_cast<npy_int32*>(pCounts+i*strideCounts)) = cnt;
            totalCount += cnt;
        }
        finishKmer();
        m_seqLen = 0;
        return bp::make_tuple(outSize,totalCount);
    }

protected:
    npy_intp m_seqLen;
    
}; // class KmerCounterPy

} // namespace MGT

using namespace boost::python;

BOOST_PYTHON_MODULE(kmersx)
{
    import_array();
    numeric::array::set_module_and_type("numpy", "ndarray");
    
    class_<MGT::KmerCounterPy,boost::noncopyable>("KmerCounter", init<int>())
        .def("process", &MGT::KmerCounterPy::process)
        .def("counts", &MGT::KmerCounterPy::counts);
}
