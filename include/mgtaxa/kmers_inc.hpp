#ifndef MGT_KMERS_INC_HPP_
#define MGT_KMERS_INC_HPP_

#ifndef MGT_KMERS_HPP_
#error THIS FILE MUST BE INCLUDED FROM kmers.hpp
#endif

#include <functional>
#include <algorithm>

namespace MGT {

////////////////////////////////////////////////////////////////
///
/// Inline implementations of functions from kmers.hpp
///
////////////////////////////////////////////////////////////////

/** Raise integer base to integer power. */

inline int ipow(int base, int power) throw (std::domain_error) {
    if (power < 0) throw std::domain_error("power must be non-negative integer");
    if (power == 0) return 1;
    if (power == 1) return base;
    if (power % 2 == 0) return ipow(base * base, power / 2);
    if (power % 2 == 1) return base * ipow(base * base, power / 2);
}


/** Convert character code to integer code. */ 

inline INuc AbcConvCharToInt::operator()(CNuc c) const {
    return m_CNucToINuc[c];
}

/** Convert character code to integer code. */

inline INuc AbcConvCharToInt::toINuc(CNuc c) const {
    return m_CNucToINuc[c];
}

/** Convert integer code to character code. */

inline CNuc AbcConvCharToInt::toCNuc(INuc i) const {
    return m_abcExt[i];
}

/** Return size of non-degenerate character alphabet*/

inline int AbcConvCharToInt::nAbc() const {
    return m_nAbc;
}

/** Return number of integer codes (with degenerate code).*/

inline int AbcConvCharToInt::nCodes() const {
    return m_nCodes;
}


/** Convert integer index of a nucleotide into integer index of a reverse-complement nucleotide.
 * @param i - integer index of a nucleotide. Must be in the valid range [0,::g_maxINuc]. 
 * Degenerate entry is converted into degenerate.*/

inline INuc AbcConvCharToInt::revCompl(INuc i) const {
    return m_iNucRevCompl[i]; 
}   

/** Print this Kmer object into ostream for debugging.*/

inline std::ostream& Kmer::print(std::ostream& out, int kmerLen) const {
    for(int i = 0; i < kmerLen; i++) {
        out << std::setw(2) << int((*this)[i]) << '\t';
    }
    return out;
}


/** Equality operator for Kmer.*/

inline bool operator==(const Kmer& x, const Kmer& y) {
    for(int i = 0; i < g_maxKmerLen; i++) {
        if( ! x[i] == y[i] ) return false;
    }
    return true;
}

inline std::ostream& operator<< (std::ostream& out, const Kmer& kmer) {
    return kmer.print(out);
}

/** Return character string representation of Kmer */

inline std::string kmerToStr(const Kmer& kmer, int kmerLen, const AbcConvCharToInt& conv) {
    std::string s(kmerLen,' ');
    for(int i = 0; i < kmerLen; i++) {
        s[i] = conv.toCNuc(kmer[i]);
    }
    return s;
}


/** Reverse-complement a k-mer.
 * @param kmerInp - input k-mer
 * @param kmerLen - length of k-mer
 * @param conv - alphabet converter
 * @param kmerOut - output k-mer, which will become the reverse-complement of kmerInp. Should not point to
 * the same location as kmerInp.
 */
 
inline void kmerToRevCompl(const Kmer& kmerInp, int kmerLen, const AbcConvCharToInt& conv, Kmer& kmerOut) {
    for(int i = 0, j = kmerLen - 1; i < kmerLen; i++, j--) {
        kmerOut[j] = conv.revCompl(kmerInp[i]);
    }
}

inline KmerState::KmerState():
    m_pData(0), 
    m_revComp(0), 
    m_isRevComp(false) {
        for(int i = 0; i < g_maxINuc; i++) {
            m_next[i] = 0;
        }
}

/** Return pointer to the next KmerState in the state machine.
 * @param pState - current state
 * @param c - nucleotide integer code
 */
  
inline PKmerState KmerStates::nextState(PKmerState pState,INuc c) {
    return pState->m_next[c];
}

/** Return pointer to the reverse complement KmerState.
 * @param pState - current state
 */

inline PKmerState KmerStates::revCompState(PKmerState pState) {
    return pState->m_revComp;
}

/** Return the state that we marked as "first" in a pair of mutually 
 * reverse-complement states.
 * @param pState points to any one of the two mutually rev-complement states.*/

inline PKmerState KmerStates::revCompStateFirst(PKmerState pState) {
    /** @todo we can avoid the branch statement by storing a pointer to
     * the "first" rev-comp state inside each KmerState object.*/
    return pState->m_isRevComp ? pState->m_revComp : pState;
}

/** Return pointer to the very first KmerState in the state machine.
 * Each k-mer count accumulation cycle starts from this state.*/ 

inline PKmerState KmerStates::firstState() {
    return &m_states[0];
}

/** Return constant pointer to the very first KmerState in the state machine.
 * Each k-mer count accumulation cycle starts from this state.*/ 


inline const PKmerState KmerStates::firstState() const {
    return const_cast<const PKmerState>(&m_states[0]);
}

/** Return total number of KmerState states in the state machine (normal and degenerate).*/ 

inline int KmerStates::numStates() const {
    return m_states.size();
}

/** Return true if this state is marked as reverse-complement.
 * This method guarantees that only one of
 * the two mutually reverse-complement states is marked as such,
 * and that it is always the same one between different invocations
 * of the program.*/
  
inline bool KmerStates::isRevComp(PKmerState pState) const {
    return pState->m_isRevComp;
}


/** Return true if a given state corresponds to a degenerate k-mer.*/

inline bool KmerStates::isDegenState(PKmerState pState) const {
    return pState - firstState() < m_blockStart.back();
}

/** Return pointer to KmerStateData object linked to the given state.
 * @param pState - current state
 * */ 

inline PKmerStateData KmerStates::getData(PKmerState pState) {
    return pState->m_pData; 
}

/** Return pointer to KmerStateData object linked to the state with a given index.
 * @param ind - index of the state
 * */ 

inline PKmerStateData KmerStates::getData(Ind ind) {
    return getData(getState(ind)); 
}

/** Link given KmerState object to the given KmerStateData object.
 * This sets only uni-directed link from KmerState to KmerStateData.
 * @param pState - pointer to state
 * @param pStateData - pointer to data
 */   

inline void KmerStates::setData(PKmerState pState,PKmerStateData pStateData) {
    pState->m_pData = pStateData;
}

/** Return reference to Kmer object for a given KmerState object.
 * @param pState - pointer to state
 * This is done in constant time based on position of pState in internal array.
 */ 

inline const Kmer& KmerStates::kmerState(const PKmerState pState) const {
    return m_kmers[pState-firstState()];
}

/** Return index of a given KmerState object that can be used to access this object later with getState().
 * @param pState - pointer to state
 * This is done in constant time based on position of pState in internal array.
 */ 

inline Ind KmerStates::indState(const PKmerState pState) const {
    return pState-firstState();
}

/** Return pointer to a  KmerState object given its index obtained earlier with indState().
 * @param ind Index of the state
 * This is done in constant time by indexing the internal array.
 */ 

inline const PKmerState KmerStates::getState(Ind ind) const {
    return const_cast<const PKmerState>(&m_states[ind]);
}

/** Return reference to a  Kmer object given its index obtained earlier with indState().
 * @param ind Index of the state
 * This is done in constant time by indexing the internal array.
 */ 

inline const Kmer& KmerStates::getKmer(Ind ind) const {
    return m_kmers[ind];
}

/** Return unique ID of the state's k-mer.
 * This is stable between invocations of the program.
 * IDs fill the half-open interval [getFirstIdState(),getLastIdState()) 
 * (for the first non-degenerate k-mer).*/

inline KmerId KmerStates::idState(PKmerState pState) const {
    return pState->m_id;
}

/** Return first state ID (for non-degenerate k-mers).*/

inline KmerId KmerStates::getFirstIdState() const {
    return m_firstIdState;
}

/** Return last used plus one state ID (for non-degenerate k-mers).*/

inline KmerId KmerStates::getLastIdState() const {
    return m_lastIdState;
}

/** Get the total number of different IDs ("number of features") - depends on the RC_POLICY.*/

inline int KmerStates::getNumIds() const {
    return getLastIdState()-getFirstIdState();
}

/** Return pointer to KmerState for a given Kmer object.
 * @param kmer - Kmer object.
 * This is done based on value of kmer, and complexity is linear 
 * in kmer length.
 */

inline PKmerState KmerStates::kmerToState(const Kmer& kmer) const {
    int ind = kmerToIndex(kmer);
    return const_cast<const PKmerState>(& m_states[ind]); 
}

/** Fill next Kmer object based on current Kmer and incoming nucleotide code.
 * If c is a degenerate code, always output the first (zero) k-mer.
 * @param currKmer - current Kmer
 * @param c - nucleotide code
 * @param nextKmer - output Kmer
 */ 

inline void KmerStates::nextKmer(const Kmer& currKmer, INuc c, Kmer& nextKmer) const {
    if( c == I_DEGEN ) {
        for(int i = 0; i < m_kmerLen; i++) {
            nextKmer[i] = I_DEGEN;
        }
    }
    else {
        int i = 0;
        for( ; i < m_kmerLen - 1; i++) {
            nextKmer[i] = currKmer[i+1];
        }
        nextKmer[i] = c;
    }
}

/** Return position of Kmer object in internal array.
 * @param kmer - Kmer object
 * This is done based on value of kmer, and complexity is linear 
 * in kmer length.
 * This is made public for the use of implementation-aware code such as KmerCounterLadder.
 */ 

inline Ind KmerStates::kmerToIndex(const Kmer& kmer) const {
    int i = 0;
    for(; i < m_kmerLen; i++) {
        if( kmer[i] != 0 ) {
            break;
        }
    }
    int ndim = m_kmerLen - i;
    int blockStart = m_blockStart[ndim];
    int blockSize = m_blockSize[ndim];
    int extent = m_pAbcConv->nAbc();
    int stride = 1;
    int ind = blockStart;
    for(; i < m_kmerLen; i++) {
        // block is indexed with non-0 codes
        ind += (kmer[i] - 1)* stride;
        stride *= extent;
    }
    return ind;
}

/** Fill the value of Kmer object from its index in the internal array.
 * @param ind - index
 * @param kmer - output Kmer object
 * This is a reverse operation to kmerToIndex() and has the same complexity.
 * @pre kmer initialized with 0s
 */

inline void KmerStates::indexToKmer(Ind ind, Kmer& kmer) const {
    int iBlock = 0;
    for(; iBlock < m_blockStart.size(); iBlock++) {
        if( m_blockStart[iBlock] > ind ) {
            break;
        }
    }
    iBlock -= 1;
    int ndim = iBlock; // iBlock == number of non-0 dimensions
    int indBlock = ind - m_blockStart[ndim];
    int extent = m_pAbcConv->nAbc();
    int iKmerNonZero = m_kmerLen - ndim;
#if ADT_DBG_LEVEL >= 4
    for( int i = 0; i < iKmerNonZero; i++ ) {
        ADT_ALWAYS(kmer[i] == I_DEGEN);
    }
#endif
    for( int i = iKmerNonZero; i < m_kmerLen; i++ ) {
        //block is indexed with non-0 codes, so we add 1:
        int coord = indBlock % extent + 1;
        kmer[i] = coord; 
        indBlock /= extent;
    }
}

/** This method is called to process the input sequence.
 * Series of calls to this method are interleaved with
 * calls to result extraction methods.
 * @param cnuc one nucleotide character value (such as 'A')*/

inline void KmerCounter::doCNuc(CNuc cnuc) {
    doINuc((*m_pAbcConv)(cnuc));
}


/** This method is called for each element of an input sequence (through
 * KmerCounter::doCNuc() wrapper method).
*/

inline void KmerCounter::doINuc(INuc inuc) {

    m_pSt = m_pStates->nextState(m_pSt,inuc);
    PKmerState pSt = m_pStates->revCompStateFirst(m_pSt);
    PKmerStateData pDat = m_pStates->getData(pSt);
    if( ! pDat ) {
        //allocate and link next data object
        pDat = &m_data[m_iDataEnd++];
        m_pStates->setData(pSt,pDat);
        pDat->setState(pSt);
    }
    pDat->count++;
}

/** Return the maximum number of unique k-mers that can be found in a sequence of a given length.
 * @param seqLen sequence length
 */

inline int KmerCounter::maxNumKmers(ULong seqLen) const {
    ULong numSeq = seqLen - m_kmerLen + 1; // an overestimate when policy is RC_MERGE
    ULong numIds = getNumIds();
    return int(numSeq<numIds ? numSeq : numIds);
}

inline int KmerCounter::numKmers() const {
    return m_iDataEnd;
}

inline ULong KmerCounter::sumKmerCounts() const {
    ULong sum = 0;
    for(int iData = 0; iData < numKmers(); iData++) {
        sum += m_data[iData].count;
    }
    return sum;
}

inline ULong KmerCounter::sumDegenKmerCounts() const {
    return m_dataDegen.count;
}

struct KmerStateLessCmp : public std::binary_function<KmerStateData,KmerStateData,bool> {
    bool operator()(const KmerStateData& x, const KmerStateData& y) {
        return x.idState() < y.idState();
    }
};

inline void KmerCounter::startKmer(bool doSort) {
    if( doSort ) {
        std::sort(m_data.begin(),m_data.begin()+numKmers(),KmerStateLessCmp());
    }
    m_iDataExtr = 0;
}

inline void KmerCounter::nextKmer() {
    KmerStateData& dat = m_data[m_iDataExtr];
    //unlink state from this data
    m_pStates->setData(dat.getState(),0);
    dat.count = 0;
    //advance to the next data point
    m_iDataExtr++;
}

inline ULong KmerCounter::getKmerCount() const {
    return m_data[m_iDataExtr].count;
}

inline int KmerCounter::getKmerId() const {
    return m_data[m_iDataExtr].idState();
    //return m_pStates->idState(m_data[m_iDataExtr].getState());
}
 
inline std::string KmerCounter::getKmerStr() const {
    return kmerToStr(m_pStates->kmerState(m_data[m_iDataExtr].getState()), m_kmerLen, *m_pAbcConv);
}

inline Ind KmerCounter::indState() const {
    return m_pStates->indState(getState());
}

inline const PKmerState KmerCounter::getState() const {
    return m_data[m_iDataExtr].getState();
}

inline void KmerCounter::finishKmer() {
    //clean up the remaining dirty KmerStateData
    for( ; m_iDataExtr < m_iDataEnd; nextKmer() )
    {}
    m_dataDegen.count = 0;
    m_iDataEnd = 0;
    m_iDataExtr = 0;
    //set current state to be the first state
    m_pSt = m_pStates->firstState();    
}

/** Return first state ID (for non-degenerate k-mers).*/

inline KmerId KmerCounter::getFirstIdState() const {
    return m_pStates->getFirstIdState();
}

/** Return last used plus one state ID (for non-degenerate k-mers).*/

inline KmerId KmerCounter::getLastIdState() const {
    return m_pStates->getLastIdState();
}

/** Get the total number of different IDs ("number of features") - depends on the RC_POLICY.*/

inline int KmerCounter::getNumIds() const {
    return m_pStates->getNumIds();
}

/** This method is called to process the input sequence.
 * Series of calls to this method are interleaved with
 * calls to result extraction methods such as counts().
 * @param cnuc one nucleotide character value (such as 'A')*/

inline void KmerCounterLadder::doCNuc(CNuc cnuc) {
    INuc inuc = (*m_pAbcConv)(cnuc);
    for(int k=m_kmerLen; k>=1; k--) {
        m_counters[k]->doINuc(inuc);
    }
}

/** Extract the counts accumulated so far and reset the internal state for the next round of accumulation.
 * @param ind output iterator for observed k-mer ID values
 * @param valObs output iterator for observed k-mer counts
 * @param valExp output iterator for expected counts for each observed k-mer
 * @param sizes output iterator for number of unique observed k-mers for each k.
 * First all unique observed k-mers for k=kmerLen will be appended to ind,valObs and valExp, 
 * and the number of written items for this k will be appended to sizes, and then it will
 * be repeated for k=kmerLen-1 and so forth, up to k=1 inclusive. Expected counts for k=1
 * will be assigned assuming equal probability. This choice is arbitrary and does not affect
 * the expected counts for higher order k-mers.*/

template<class IterInd,class IterVal>
void KmerCounterLadder::counts(IterVal valObs,IterVal valExp,IterInd ind,IterInd sizes) {
    //We assign expected values for 1-mers as 1/4*length. The user is free to
    //override these post-factum as nothing that is computed here depends on them.
    //
    //Expected count of k-mer: N_expect(k) = N_observed(k-1)*p(1)
    KmerCounter& one = *m_counters[1];
    KmerStates& oneStates = one.getStates();
    ULong sumOne = one.sumKmerCounts();
    // k==1 should be treated differently - we achieve it by
    // inserting condition within the inner loop to avoid
    // duplicating the code
    for(int k=m_kmerLen; k>0; k--) {
        KmerCounter& top = *m_counters[k];
        // bot should not be used if k==1
        KmerCounter& bot = *m_counters[k>1?k-1:k];
        KmerStates& topStates = top.getStates();
        KmerStates& botStates = bot.getStates();
        // should not be used if k==1
        indvec& topBotDep = m_topBotDep[k];
        indvec& topOneDep = m_topOneDep[k];
        int n = top.numKmers();
        top.startKmer();
        for(int i = 0; i < n; i++,top.nextKmer()) {
            double countExpTop = 0;
            if(k>1) {
                Ind indStateTop = top.indState();
                Ind indStateBot = topBotDep[indStateTop];
                Ind indStateOne = topOneDep[indStateTop];
                PKmerStateData pStateDataBot = botStates.getData(indStateBot);
                PKmerStateData pStateDataOne = oneStates.getData(indStateOne);
                //if we have seen a k-mer (top), we should have seen all prefixes (bottom),
                //and so bottom state object has its data payload linked
                ADT_ALWAYS(pStateDataBot);
                ADT_ALWAYS(pStateDataOne);
                countExpTop = pStateDataBot->count * (double(pStateDataOne->count)/sumOne);
            }
            else {
                countExpTop = double(sumOne)/n;
            }
            *ind = top.getKmerId(); ++ind;
            *valObs = top.getKmerCount(); ++valObs;
            *valExp = countExpTop; ++valExp;
        }
        top.finishKmer();
        *sizes = n; ++sizes;
    }
}

/** Store into the output iterator the number of unique k-mers found for each k and return their sum.
 * @param output iterator - will hold a number of unique k-mer for k,k-1,...,1
 */

template<class IterInd>
int KmerCounterLadder::numKmers(IterInd sizes) const {
    int n_sum = 0;
    for(int k=m_kmerLen; k>0; k--) {
        int n = m_counters[k]->numKmers();
        n_sum += n;
        *sizes = n; ++sizes;
    }
    return n_sum;
}

/** Store into the output iterator the maximum number of unique k-mers that can be found for each k for a given sequence length.
 * @param seqLen sequence length
 * @param output iterator - will hold a number of unique k-mer for k,k-1,...,1
 * @return sum of values stored into sizes
 */

template<class IterInd>
int KmerCounterLadder::maxNumKmers(ULong seqLen, IterInd sizes) const {
    int n_sum = 0;
    for(int k=m_kmerLen; k>0; k--) {
        int n = m_counters[k]->maxNumKmers(seqLen);
        n_sum += n;
        *sizes = n; ++sizes;
    }
    return n_sum;
}

} // namespace MGT

#endif // MGT_KMERS_INC_HPP_
