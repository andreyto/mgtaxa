#include "mgtaxa/kmers.hpp"

#define PY_ARRAY_UNIQUE_SYMBOL PyArrayHandle

#include <boost/python.hpp>

#include "mgtaxa/py_num_util.hpp"

#include <cstdio>
#include <zlib.h>

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

    protected:
        enum RES_TYPE { RES_COUNTS, RES_BITS, RES_FREQUENCES };

public:
	
    KmerCounterPy(int kmerLen,KmerId firstIdState=1,bool doSort=true) :
	KmerCounter(kmerLen,0,firstIdState),
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
        return tp_counts<npy_int32>(counts,indices,NPY_INT32,RES_COUNTS);
    }

    bp::tuple bits(bpn::array counts, bpn::array indices) {
        return tp_counts<npy_int32>(counts,indices,NPY_INT32,RES_BITS);
    }

    bp::tuple frequences(bpn::array counts, bpn::array indices) {
        return tp_counts<npy_float32>(counts,indices,NPY_FLOAT32,RES_FREQUENCES);
    }

    void setSort(bool doSort) {
        m_doSort = doSort;
    }

    bp::tuple getRangeIdState() const {
        return bp::make_tuple(getFirstIdState(),getLastIdState());
    }

protected:
    
    template<typename CTypeCounts>
    bp::tuple tp_counts(bpn::array counts, 
            bpn::array indices,
            NPY_TYPES npyTypeCounts,
            RES_TYPE resType) {
        bool divide = false;
        if( resType == RES_FREQUENCES ) {
            divide = true;
        }
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
            if( resType == RES_BITS ) {
                cnt = cnt > 0? 1 : 0;
            }
            int id = getKmerId();
            (*reinterpret_cast<npy_int64*>(pIndices+i*strideIndices)) = id;
            (*reinterpret_cast<CTypeCounts*>(pCounts+i*strideCounts)) = cnt;
            totalCount += cnt;
        }
        finishKmer();
        m_seqLen = 0;
        if( divide ) {
            CTypeCounts divider = CTypeCounts(totalCount)/100;
            for(npy_intp i = 0; i < outSize; i++) {
                (*reinterpret_cast<CTypeCounts*>(pCounts+i*strideCounts)) /= divider;
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
	
    KmerCounterWriterPy(int kmerLen,KmerId firstIdState=1,bool doSort=true) :
	KmerCounterPy(kmerLen,firstIdState,doSort),
    m_file(0)
	{
        //ADT_LOG << ADT_OUTVAR(kmerLen) << '\n';
    }

    void openOutput(const std::string& fileName) {
        if( m_file ) {
            std::fclose(m_file);
            //gzclose(m_file);
        }
        m_file = std::fopen(fileName.c_str(),"w");
        //m_file = gzopen(fileName.c_str(),"w");
    }

    void closeOutput() {
        if( m_file ) {
            std::fclose(m_file);
            //gzclose(m_file);
            m_file = 0;
        }
    }


    ULong frequencesWriteSvmSparseTxt(int label,int maxDegen) {
        ULong outSize = numKmers();
        //double sumCounts = double(sumKmerCounts())/100;
        ULong nDegen = sumDegenKmerCounts();
        ULong nWritten = 0; 
        startKmer(true);
        if( nDegen < maxDegen ) {
            std::fprintf(m_file,"%i",label);
            //gzprintf(m_file,"%i",label);
            for(ULong i = 0; i < outSize; i++,nextKmer()) {
                std::fprintf(m_file," %i:%i",getKmerId(),getKmerCount() > 0? 1 : 0);
                //std::fprintf(m_file," %i:%f",getKmerId(),getKmerCount()/sumCounts);
                //gzprintf(m_file," %i:%f",getKmerId(),getKmerCount()/sumCounts);
            }
            std::fprintf(m_file,"\n");
            //gzprintf(m_file,"\n");
            nWritten++;
        }
        finishKmer();
        return nWritten;
    }

    virtual ~KmerCounterWriterPy() {
        closeOutput();
    }

protected:
    std::FILE * m_file;
    //gzFile m_file;

}; // class KmerCounterWriterPy

/*
 
class SvmSparseFeatureWriterTxt:
    
    def __init__(self,out):
        if not hasattr(out,'write'):
            out = open(out,'w', buffering=1024*1024)
        self.out = out
        self.nOut = 0
        
    def close(self):
        self.out.close()

    def write(self,label,values,indices):
        self.out.write("%d " % (label,))
        strItems = ' '.join( ( "%d:%f" % item for item in izip(indices,values) ) )
        self.out.write(strItems)
        self.out.write("\n")
        self.nOut += 1

    def numRec(self):
        return self.nOut
 */

class SvmSparseFeatureWriterTxt : boost::noncopyable {
    public:
    SvmSparseFeatureWriterTxt(const std::string& fileName) {
        m_file = std::fopen(fileName.c_str(),"w");
        m_nOut = 0;
    }
    ~SvmSparseFeatureWriterTxt() {
        this->close();
    }
    void close() {
        if( m_file ) {
            std::fclose(m_file);
            m_file = 0;
        }
    }

    template<typename CTypeValues>
    void tp_writePy(NPY_TYPES npyTypeValues,
                    int label, 
                    const bpn::array& values, 
                    const bpn::array& indices) {
        
        npy_intp sizeValues = 0, sizeIndices = 0, strideValues = 0, strideIndices = 0;
        char * pValues = getCData1D<CTypeValues>(values,npyTypeValues,sizeValues,strideValues);
        char * pIndices = getCData1D<npy_int64>(indices,NPY_INT64,sizeIndices,strideIndices);

        if( sizeValues != sizeIndices){
            PyErr_SetString(PyExc_ValueError, "Expected value and index arrays of equal length");
            bp::throw_error_already_set();
        }

        std::fprintf(m_file,"%i",label);
        for(ULong i = 0; i < sizeValues; i++) {
            npy_int64 id = (*reinterpret_cast<npy_int64*>(pIndices+i*strideIndices));
            double v = double(*reinterpret_cast<CTypeValues*>(pValues+i*strideValues));
            std::fprintf(m_file," %i:%g",id,v);
        }
        std::fprintf(m_file,"\n");

        m_nOut++;
    }

    void writePy(int label, 
                 const bpn::array& values, 
                 const bpn::array& indices) {

        NPY_TYPES npyTypeValues = num_util::type(values);

        switch(npyTypeValues) {
            case NPY_FLOAT32:
                tp_writePy<npy_float32>(npyTypeValues,label,values,indices);
                break;
            case NPY_FLOAT64:
                tp_writePy<npy_float64>(npyTypeValues,label,values,indices);
                break;
            case NPY_INT32:
                tp_writePy<npy_int32>(npyTypeValues,label,values,indices);
                break;
            case NPY_INT64:
                tp_writePy<npy_int64>(npyTypeValues,label,values,indices);
                break;
            case NPY_INT16:
                tp_writePy<npy_int16>(npyTypeValues,label,values,indices);
                break;
            case NPY_BOOL:
                tp_writePy<npy_bool>(npyTypeValues,label,values,indices);
                break;
            default:
                PyErr_SetString(PyExc_ValueError, "Unsupported type for the values array");
                bp::throw_error_already_set();
        }
    }

    int numRec() const {
        return m_nOut;
    }

    protected:
    std::FILE * m_file;
    int m_nOut;
};

} // namespace MGT

using namespace boost::python;

BOOST_PYTHON_MODULE(kmersx)
{
    import_array();
    numeric::array::set_module_and_type("numpy", "ndarray");
    
    class_<MGT::KmerCounterWriterPy,boost::noncopyable>("KmerCounter", 
            init<int,optional<MGT::KmerId,bool> >())
        .def("process", &MGT::KmerCounterWriterPy::process)
        .def("counts", &MGT::KmerCounterWriterPy::counts)
        .def("numKmers", &MGT::KmerCounterWriterPy::numKmers)
        .def("sumDegenKmerCounts", &MGT::KmerCounterWriterPy::sumDegenKmerCounts)
        .def("frequences", &MGT::KmerCounterWriterPy::frequences)
        .def("bits", &MGT::KmerCounterWriterPy::bits)
        .def("getFirstIdState", &MGT::KmerCounterWriterPy::getFirstIdState)
        .def("getLastIdState", &MGT::KmerCounterWriterPy::getLastIdState)
        .def("getRangeIdState", &MGT::KmerCounterWriterPy::getRangeIdState)
        .def("setSort",&MGT::KmerCounterWriterPy::setSort)
        .def("openOutput", &MGT::KmerCounterWriterPy::openOutput)
        .def("closeOutput", &MGT::KmerCounterWriterPy::closeOutput)
        .def("frequencesWriteSvmSparseTxt", &MGT::KmerCounterWriterPy::frequencesWriteSvmSparseTxt);
    
    class_<MGT::SvmSparseFeatureWriterTxt,boost::noncopyable>("SvmSparseFeatureWriterTxt", 
            init<const std::string&>())
        .def("write", &MGT::SvmSparseFeatureWriterTxt::writePy)
        .def("close", &MGT::SvmSparseFeatureWriterTxt::close)
        .def("numRec", &MGT::SvmSparseFeatureWriterTxt::numRec);
}

