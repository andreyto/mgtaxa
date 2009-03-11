#include "mgtaxa/kmers.hpp"

#define PY_ARRAY_UNIQUE_SYMBOL PyArrayHandle

#include <boost/python.hpp>
#include <boost/iterator/iterator_facade.hpp>

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

/** STL compliant random access iterator for 1D Numpy arrays exposed as boost:python::array.
 * This iterates using stride and data pointer and thus more efficient than
 * the generic ND implementation with PyArrayIterObject.
 * It will correctly handle arrays which have stride not in multiples of data 
 * element size, such as column views of record arrays.*/

template<typename TVal>
class npy_iterator_1d
  : public boost::iterator_facade<
    npy_iterator_1d<TVal>
    , TVal
    , boost::random_access_traversal_tag
    >
{
    public:
    typedef npy_iterator_1d<TVal> Self;

    npy_iterator_1d()
    : m_data(0), m_stride(0) {}

    /** Constructor.
     * @param data pointer to the start of the current array element
     * @param stride distance between the starts of consequitive elements in bytes, as returned by PyArray_STRIDE(a,0)
     */
    explicit npy_iterator_1d(char* data,npy_intp stride)
    : m_data(data), m_stride(stride) {}

    private:
    friend class boost::iterator_core_access;

    void increment() { m_data += m_stride; }
    
    void decrement() { m_data -= m_stride; }
    
    void advance(typename Self::difference_type n) { m_data += n*m_stride; }

    typename Self::difference_type distance_to(const Self& other) const 
    { 
        return (other->m_data-this->m_data)/this->m_stride;
    }

    bool equal( const Self& other) const
    {
       return this->m_data == other.m_data;
    }

    TVal& dereference() const { return (*reinterpret_cast<TVal*>(m_data)); }

    char* m_data;
    npy_intp m_stride;
};


/** Wrapper around boost::python::array Numpy 1D array object that can generate STL compliant 1D iterators.*/

template<typename TVal>
class npy_array_1d_wrapper {
    
    public:
    typedef npy_iterator_1d<TVal> iterator_type;
    
    npy_array_1d_wrapper(bpn::array a):
        m_a(a) 
    {
        num_util::check_rank(m_a,1);
        //num_util::check_type(a,npyType);
        /**@todo add type compatibility check - will require creating a set of template trait specializations*/
        if( sizeof(TVal) != PyArray_ITEMSIZE(m_a.ptr()) ){
            PyErr_SetString(PyExc_ValueError, "C type size does not match NumPy array item size");
            bp::throw_error_already_set();
        }
    }

    iterator_type begin() {
        return iterator_type(reinterpret_cast<char*>(num_util::data(m_a)),PyArray_STRIDE(m_a.ptr(),0));
    }
    
    iterator_type end() {
        return begin()+num_util::size(m_a);
    }

    bpn::array get_array() {
        return m_a;
    }

    npy_intp size() const {
        return num_util::size(m_a);
    }

    protected:
    bpn::array m_a;
};


class KmerCounterPy : public KmerCounter {

    protected:
        enum RES_TYPE { RES_COUNTS, RES_BITS, RES_FREQUENCES };

public:
	
    KmerCounterPy(int kmerLen,RC_POLICY revCompPolicy=RC_MERGE,KmerId firstIdState=1,bool doSort=true) :
	KmerCounter(kmerLen,0,revCompPolicy,firstIdState),
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
            PyErr_SetString(PyExc_ValueError, "Output size exceeds array bounds");
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
            CTypeCounts divider = CTypeCounts(totalCount);
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
	
    KmerCounterWriterPy(int kmerLen,RC_POLICY revCompPolicy=RC_MERGE,KmerId firstIdState=1,bool doSort=true) :
	KmerCounterPy(kmerLen,revCompPolicy,firstIdState,doSort),
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


class KmerCounterLadderPy : public KmerCounterLadder {

public:
	
    KmerCounterLadderPy(int kmerLen,RC_POLICY revCompPolicy=RC_MERGE) :
	KmerCounterLadder(kmerLen,0,revCompPolicy)
	{
        //ADT_LOG << ADT_OUTVAR(kmerLen) << '\n';
    }
	
	void processPy(bpn::array seq) {
        npy_intp nSeq = 0, stride = 0;
        const char *pSeq = getCData1D<char>(seq,NPY_CHAR,nSeq,stride);
        //ADT_LOG <<ADT_OUTVAR(nSeq) << ADT_OUTVAR(pSeq) << ADT_OUTVAR(stride) << '\n';
        for(npy_intp i = 0; i < nSeq; i++) {
            doCNuc(pSeq[i*stride]);
        }
	}

    void countsPy(bpn::array valObs,bpn::array valExp,bpn::array ind,bpn::array sizes) {
        typedef npy_array_1d_wrapper<npy_float32> npyw_val;
        typedef npy_array_1d_wrapper<npy_int64> npyw_ind;
        npyw_val valObsW(valObs), valExpW(valExp);
        npyw_ind indW(ind), sizesW(sizes);
        ADT_ALWAYS(sizesW.size()>=numSubFeatures());
        int n_feat = numKmers(sizesW.begin());
        ADT_ALWAYS(valObsW.size()>=n_feat && valExpW.size()>=n_feat && indW.size()>=n_feat);
        counts(valObsW.begin(),valExpW.begin(),indW.begin(),sizesW.begin());
    }

    bpn::array maxNumKmersPy(ULong seqLen) const {
        typedef npy_array_1d_wrapper<npy_int64> npyw_ind;
        bpn::array sizes = num_util::makeNum(m_kmerLen,NPY_INT64);
        maxNumKmers(seqLen,npyw_ind(sizes).begin());
        return sizes;
    }

};

} // namespace MGT

using namespace boost::python;

BOOST_PYTHON_MODULE(kmersx)
{
    import_array();
    numeric::array::set_module_and_type("numpy", "ndarray");
    
    class_<MGT::KmerCounterWriterPy,boost::noncopyable>("KmerCounter", 
            init<int,optional<MGT::RC_POLICY,MGT::KmerId,bool> >())
        .def("process", &MGT::KmerCounterWriterPy::process)
        .def("counts", &MGT::KmerCounterWriterPy::counts)
        .def("maxNumKmers", &MGT::KmerCounterWriterPy::maxNumKmers)
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
    
    enum_<MGT::RC_POLICY>("RC_POLICY")
        .value("MERGE", MGT::RC_MERGE)
        .value("DIRECT", MGT::RC_DIRECT)
        ;

    class_<MGT::KmerCounterLadderPy,boost::noncopyable>("KmerCounterLadder", 
            init<int,optional<MGT::RC_POLICY> >())
        .def("process", &MGT::KmerCounterLadderPy::processPy)
        .def("counts", &MGT::KmerCounterLadderPy::countsPy)
        .def("maxNumKmers", &MGT::KmerCounterLadderPy::maxNumKmersPy);
}

