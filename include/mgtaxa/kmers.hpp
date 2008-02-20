#ifndef MGT_KMERS_HPP_
#define MGT_KMERS_HPP_

#include "mgtaxa/types.hpp"
#include "mgtaxa/exceptions.hpp"

#include <vector>
#include <iostream>
#include <cmath>
#include <string>

namespace MGT {

/** Alphabet convertor from one-letter character to integer index.*/
 
class AbcConvCharToInt {

	public:

	AbcConvCharToInt(const std::string& abc);
	~AbcConvCharToInt();
	
	INuc operator()(CNuc c) const;
	CNuc toCNuc(INuc i) const;
	
	
	protected:

	/**If abc is 'ACTG', then m_CNucToINuc['A'] => 1, m_CNucToINuc['C'] => 2 and so on. 0 is reserved for 'N'*/
	const CNuc * m_CNucToINuc;
	
	/**Extended alphabet, e.g. - 'NACTG'*/ 
	std::string m_abcExt;
	
	/**number of "real" alphabet symbols*/
	int m_nAbc;
	/**total number of encoded alphabet symbols (m_nAbc + 1)*/
	int m_nINuc;

};		


/** Class that describes a k-mer in integer index representation.
 * We expect to create many k-mer objects for the same k, so k
 * is stored not inside this class but somewhere else (in the container).*/ 

class Kmer {
	public:
	Kmer() {
		for(int i = 0; i < g_maxKmerLen; i++) { m_data[i] = 0; }
	}
	INuc m_data[g_maxKmerLen];
	inline INuc& operator[](int i) { return m_data[i]; }
	inline const INuc& operator[](int i) const { return m_data[i]; }	
};


std::string kmerToStr(const Kmer& kmer, int kmerLen, const AbcConvCharToInt& conv);  

/** Forward declaration of payload class for  KmerState.*/

class KmerStateData;

typedef KmerStateData *PKmerStateData;

/** Type for unique, stable between program invocations ID of a k-mer.
 * For dense  KmerStates implementation, this is defined as int.
 * For sparce (e.g. hash) representations, it might be k-mer itself.*/
 
typedef int KmerId;

/** Class that holds a counter for one k-mer along with jump pointers for next k-mers.*/

class KmerState {
	public:

	KmerState();

	typedef KmerState* PKmerState;

	public:

	/** Array of "edges" - a pointer to another  KmerState for each extended alphabet symbol.
	 * Declared public only for access by KmerStates class.*/
	PKmerState m_next[g_maxINuc];

	/** Pointer to payload object.
	 * Declared public only for access by  KmerStates class.*/
	PKmerStateData m_pData;
	
	/** Pointer to a reverse-complement  KmerState object.
	 * Declared public only for access by  KmerStates class.*/
	PKmerState m_revComp;	

	/** true if we marked this state as reverse-complement (for KmerStates::isRevComp()).
	 * Declared public only for access by  KmerStates class.*/
	bool m_isRevComp;
	
	/** ID of corresponding k-mer (for KmerStates::idState()).
	 * Declared public only for access by  KmerStates class.*/ 
	KmerId m_id;	
		
};

/** Pointer type for  KmerCounter */

typedef KmerState::PKmerState PKmerState;


/** Payload data class for  KmerState.
 * Separate class for payload allows lazy allocation of that data.
 * The main reason is to have as few payload objects as we actually saw k-mers
 * in the input. This in turn will allow extraction of k-mer counts in time
 * which does not depend on the size of k-mer space (k-mer length).*/

class KmerStateData {
	public:

	KmerStateData():
	m_pState(0), 
	count(0)
	{}
	
	PKmerState getState() {
		return m_pState;
	}
	
	void setState(PKmerState pState) {
		m_pState = pState;
	}

	/** k-mer count.*/
	ULong count;
	
	protected:
	
	/** Pointer to the owning  KmerState object.*/
	PKmerState m_pState;
	
};


/** State machine for k-mers along with a corresponding payload object
 * implemented in  KmerStates and  KmerCounter classes respectively. 
 * The current implementation is optimized for the following use case:
 * <ol>
 * <li> relatively short k-mers (k < 10)</li>
 * <li> @anchor AnchShortSeqInp k-mer counts will be computed for many short chunks of sequence, which assumes
 * the following set of operations</li>
 * <ul> <li>reset counts</li> <li>stream through sequence chunk</li> <li>report the counts</li> 
 * <li>repeat</li></ul>
 * Typical length of sequence chunk (~1000 nuc) is less than the number of
 * possible k-mers (~16,000 for k = 7).
 * <li> k-mer counts will be extracted only for non-redundant subset of 
 * reverse complement k-mers (~8,000 for k = 7)</li>
 * </ol>
 * In view of @ref AnchShortSeqInp "2", we will keep the "payload" for each state node separately,
 * referenced by the pointer to an instance of KmerStateData class. The payload object will 
 * contain the back reference
 * to its state and the count, and the payloads can be allocated as needed, from a pre-allocated
 * array. The size of payload array is bound by the max chunk length x 2 (direct plus reverse-complement
 * sequence lengths), or the total number of possible k-mers for a given k, whichever is less.
 * Therefore, we pre-allocate storage for the payload array of with the total number of possible k-mers.
 * 
 * The current implementation uses dense k-mer state machine (all possible k-mers are
 * preallocated and states are initialized) with sparse payload: all states have initially 
 * their payloads set to NULL. As we stream through the input sequence, we move from the current
 * state to the next one. At each step, we take the first of the reverse complement pair of states,
 * and check if its payload is NULL. If it is, we allocate the new payload object. Then, we 
 * increment the count in the payload object.
 *
 * When counts are to be extracted, we have to loop only through the payload objects for those k-mers 
 * that were actually seen and represent the non-redundant subset of all reverse-complement pairs.  
 * 
 * After construction, this implementation guarantees the following complexity during
 * the sequence processing: O(1) on k-mer length; O(min(L,N)) where L is the actual sequence chunk
 * length and N is the number of possible k-mers.*/ 

class KmerStates {
	protected:
	typedef std::vector<KmerState> StateArray;
	typedef std::vector<Kmer> KmerArray; 
	public:
		KmerStates(int kmerLen,int nAbc);
		PKmerState nextState(PKmerState pState,INuc c);
		PKmerState revCompState(PKmerState pState);
		PKmerState revCompStateFirst(PKmerState pState);		
		bool isRevComp(PKmerState pState);
		PKmerState firstState();
		int numStates() const;
		KmerId idState(PKmerState pState);
		PKmerStateData getData(PKmerState pState);
		void setData(PKmerState pState,PKmerStateData pStateData);
		const Kmer& kmerState(PKmerState pState) const;
	protected:
		PKmerState kmerToState(const Kmer& kmer); 
		void nextKmer(const Kmer& currKmer, INuc c, Kmer& nextKmer);
		int kmerToIndex(const Kmer& kmer);
		void indexToKmer(int ind,Kmer& kmer);
		void initAllKmers();

	protected:
	/**number of "real" alphabet symbols.*/
	int m_nAbc;
	/**total number of encoded alphabet symbols - (m_nAbc + 1)*/
	int m_nCodes;
	int m_kmerLen;
	std::vector<int> m_blockSize;
	std::vector<int> m_blockStart;
	
	StateArray m_states;
	KmerArray m_kmers;
	
};


/** Class that computes k-mer counts of incoming nucleotide sequence.
It processes the input sequence and extracts the counts.
Internally, it maintains the KmerStateData payload data and moves through
states of the KmerStates state machine in response to incoming nucleotides.*/
 
class KmerCounter {
	
	public:
	
	KmerCounter(int kmerLen, const std::string& abc);
	~KmerCounter();
	
	void doCNuc(CNuc cnuc);
	void doINuc(INuc inuc);

	/** @name Interface to extract the results.
 	* The following set of methods defines the result extraction protocol.
 	* They must be called in a strict order because they change the internal state
 	* of KmerCounter object. The interface is designed for efficiency and flexibility
 	* of the caller.
 	* After doCNuc() has been called any number of times, the extraction is done as in the
 	* following code sample, after which doCNuc() can be called again to accumulate new counts.
 	* @code
 	* int n = o.numKmers();
 	* o.startKmer();
 	* for(int i = n; i < n; i++,o.nextKmer()) {
 	*     cout << o.getKmerId() << ":" << o.getKmerCount() << "\n";
 	* }
 	* o.finishKmer();
 	* @endcode
 	* nextKmer() can be called less than numKmers() times.
 	* */
	/*@{*/ 	
	/** Return number of k-mers found so far.*/
	int numKmers() const;
	/** Prepare internal state for result extraction.*/
	void startKmer();
	/** Advance internal state to extract next k-mer results.*/
	void nextKmer();
	/** Accessor to get count value from the currently extracted k-mer.*/
	ULong getKmerCount() const;
	/** Accessor to get Id from the currently extracted k-mer.*/
	int getKmerId() const;
	/** Accessor to get k-mer string (such as 'ACCCT') for the currently extracted k-mer.*/ 
	std::string getKmerStr() const;
	/** Finalize result extraction cycle - new series of doCNuc() calls can be done afterwards.*/
	void finishKmer();
	/*@}*/ 
	
	protected:

	AbcConvCharToInt * m_pAbcConv;
	
	PKmerState m_pSt;
	KmerStates *m_pStates;
	
	int m_kmerLen;
	
	/** Preallocated array of  KmerStateData objects.
	 * Array size is equal to the total number of states.
	 * This guarantees that there is always enough data elements
	 * regardless of the input sequence length.*/
	std::vector<KmerStateData> m_data;
	
	/** Index of the first unused element in ::m_data.*/
	int m_iDataEnd;
	
	/** Index of m_data element that is currently being extracted.*/
	int m_iDataExtr; 
};


////////////////////////////////////////////////////////////////
///
/// Inline implementations
///
////////////////////////////////////////////////////////////////

inline INuc AbcConvCharToInt::operator()(CNuc c) const {
	return m_CNucToINuc[c];
}

inline CNuc AbcConvCharToInt::toCNuc(INuc i) const {
	return m_abcExt[i];
}

inline std::string kmerToStr(const Kmer& kmer, int kmerLen, const AbcConvCharToInt& conv) {
	std::string s(kmerLen);
	for(int i = 0; i < kmerLen; i++) {
		s[i] = conv.toCNuc(kmer[i]);
	}
	return s;
}

inline KmerState::KmerState():
	m_pData(0), 
	m_revComp(0), 
	m_isRevComp(false) {
		for(int i = 0; i < g_maxINuc; i++) {
			m_next[i] = 0;
		}
}

inline PKmerState KmerStates::nextState(PKmerState pState,INuc c) {
	return pState->m_next[c];
}

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

inline PKmerState KmerStates::firstState() {
	return &m_states[0];
}

inline int KmerStates::numStates() const {
	return m_states.size();
}

/** Return true if this state is marked as reverse-complement.
 * This method guarantees that only one of
 * the two mutually reverse-complement states is marked as such,
 * and that it is always the same one between different invocations
 * of the program.*/
  
inline bool KmerStates::isRevComp(PKmerState pState) {
	return pState->m_isRevComp;
}

inline PKmerStateData KmerStates::getData(PKmerState pState) {
	return pState->m_pData; 
}

inline void KmerStates::setData(PKmerState pState,PKmerStateData pStateData) {
	pState->m_pData = pStateData;
}

inline const Kmer& KmerStates::kmerState(PKmerState pState) const {
	return m_kmers[pState-firstState()];
}

/** Return unique ID of the state's k-mer.
 * This is stable between invocations of the program.*/

inline KmerId KmerStates::idState(PKmerState pState) {
	return pState->m_id;
}

inline PKmerState KmerStates::kmerToState(const Kmer& kmer) {
	int ind = kmerToIndex(kmer);
	return & m_states[ind]; 
}

inline int KmerStates::kmerToIndex(const Kmer& kmer) {
	int i = 0;
	for(; i < m_kmerLen; i++) {
		if( kmer[i] != 0 ) {
			break;
		}
	}
	int ndim = m_kmerLen - i;
	int blockStart = m_blockStart[ndim];
	int blockSize = m_blockSize[ndim];
	int extent = m_nAbc;
	int stride = 1;
	int ind = blockStart;
	for(; i < m_kmerLen; i++) {
		ind += kmer[i] * stride;
		stride *= extent;
	}
	return ind;
}


inline void KmerStates::indexToKmer(int ind, Kmer& kmer) {
	int iBlock = 0;
	for(; iBlock < m_kmerLen; iBlock++) {
		if( m_blockStart[iBlock] > ind ) {
			break;
		}
	}
	iBlock -= 1;
	int ndim = iBlock; // iBlock == number of non-0 dimensions
	int indBlock = ind - m_blockStart[ndim];
	int extent = m_nAbc;
	int stride = m_nAbc; 
	for( int i = m_kmerLen - ndim; i < m_kmerLen; i++ ) {
		kmer[i] = indBlock % stride;
		stride *= extent;
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
 * @todo make sure that reverse-complement pointers work correctly
 * for degenerate k-mer states and palyndromic states.
 * @todo for degenerates, we can create a separate dummy KmerStateData
 * object, and point all such states to it. That will 
 * leave only valid output counts in the main data array.*/

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


inline int KmerCounter::numKmers() const {
	return m_iDataEnd;
}

inline void KmerCounter::startKmer() {
	m_iDataExtr = 0;
}

inline void KmerCounter::nextKmer() {
	KmerStateData& dat = m_data[m_iDataExtr];
	//unlink state from this data
	m_pStates->setData(dat.getState(),0);
	//this is not actually needed, but just in case,
	//unlink data from state
	dat.setState(0);
	dat.count = 0;
	//advance to the next data point
	m_iDataExtr++;
}

inline ULong KmerCounter::getKmerCount() const {
	return m_data[m_iDataExtr].count;
}

inline int KmerCounter::getKmerId() const {
	return m_pStates->idState(m_data[m_iDataExtr].getState())
}
 
inline std::string KmerCounter::getKmerStr() const {
	return kmerToStr(m_pStates->kmerState(m_data[m_iDataExtr].getState()), m_kmerLen, *m_pAbcConv);
}

inline void KmerCounter::finishKmer() {
	//clean up the remaining dirty KmerStateData
	for( ; m_iDataExtr < m_iDataEnd; nextKmer() )
	{}
	m_iDataEnd = 0;
	m_iDataExtr = 0;
}


} // namespace MGT

#endif /*MGT_KMERS_HPP_*/
