#include "mgtaxa/kmers.hpp"

#define PY_ARRAY_UNIQUE_SYMBOL PyArrayHandle

#include <boost/python.hpp>

#include "mgtaxa/py_num_util.hpp"

#include <cstdio>

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
	
    KmerCounterPy(int kmerLen,bool doSort) :
	KmerCounter(kmerLen),
    m_seqLen(0),
    m_doSort(doSort)
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
        return tp_counts<npy_int32>(counts,indices,NPY_INT32,false);
    }

    bp::tuple frequences(bpn::array counts, bpn::array indices) {
        return tp_counts<npy_float32>(counts,indices,NPY_FLOAT32,true);
    }

    void setSort(bool doSort) {
        m_doSort = doSort;
    }

protected:
    
    template<typename CTypeCounts>
    bp::tuple tp_counts(bpn::array counts, 
            bpn::array indices,
            NPY_TYPES npyTypeCounts,
            bool divide) {
        //ADT_LOG << ADT_OUTVAR(m_seqLen) << '\n';
        npy_intp totalCount = 0;
        npy_intp sizeCounts = 0, sizeIndices = 0, strideCounts = 0, strideIndices = 0;
        char * pCounts = getCData1D<CTypeCounts>(counts,npyTypeCounts,sizeCounts,strideCounts);
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
        startKmer(m_doSort);
        for(npy_intp i = 0; i < outSize; i++,nextKmer()) {
            npy_intp cnt = getKmerCount();
            int id = getKmerId();
            (*reinterpret_cast<npy_int64*>(pIndices+i*strideIndices)) = id;
            (*reinterpret_cast<CTypeCounts*>(pCounts+i*strideCounts)) = cnt;
            totalCount += cnt;
        }
        finishKmer();
        m_seqLen = 0;
        if( divide ) {
            for(npy_intp i = 0; i < outSize; i++) {
                (*reinterpret_cast<CTypeCounts*>(pCounts+i*strideCounts)) /= totalCount;
            }
        }
        return bp::make_tuple(outSize,totalCount);
    }

protected:
    npy_intp m_seqLen;
    bool m_doSort;
}; // class KmerCounterPy


class KmerCounterWriterPy : public KmerCounterPy {

public:
	
    KmerCounterWriterPy(int kmerLen,bool doSort) :
	KmerCounterPy(kmerLen,doSort),
    m_file(0)
	{
        //ADT_LOG << ADT_OUTVAR(kmerLen) << '\n';
    }

    void openOutput(const std::string& fileName) {
        if( m_file ) {
            std::fclose(m_file);
        }
        m_file = std::fopen(fileName.c_str(),"w");
    }

    void closeOutput() {
        if( m_file ) {
            std::fclose(m_file);
            m_file = 0;
        }
    }


    void frequencesWriteSvmSparseTxt(int label) {
        std::fprintf(m_file,"%i",label);
        ULong outSize = numKmers();
        double sumCounts = double(sumKmerCounts()); 
        startKmer(true);
        for(ULong i = 0; i < outSize; i++,nextKmer()) {
            std::fprintf(m_file," %i:%f",getKmerId(),getKmerCount()/sumCounts);
        }
        finishKmer();
        std::fprintf(m_file,"\n");
    }

    virtual ~KmerCounterWriterPy() {
        closeOutput();
    }

protected:
    std::FILE * m_file;

}; // class KmerCounterWriterPy

} // namespace MGT

using namespace boost::python;

BOOST_PYTHON_MODULE(kmersx)
{
    import_array();
    numeric::array::set_module_and_type("numpy", "ndarray");
    
    class_<MGT::KmerCounterWriterPy,boost::noncopyable>("KmerCounter", init<int,bool>())
        .def("process", &MGT::KmerCounterWriterPy::process)
        .def("counts", &MGT::KmerCounterWriterPy::counts)
        .def("frequences", &MGT::KmerCounterWriterPy::frequences)
        .def("setSort",&MGT::KmerCounterWriterPy::setSort)
        .def("openOutput", &MGT::KmerCounterWriterPy::openOutput)
        .def("closeOutput", &MGT::KmerCounterWriterPy::closeOutput)
        .def("frequencesWriteSvmSparseTxt", &MGT::KmerCounterWriterPy::frequencesWriteSvmSparseTxt);
}
