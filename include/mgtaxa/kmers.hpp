#ifndef MGT_KMERS_HPP_
#define MGT_KMERS_HPP_

#include "mgtaxa/types.hpp"
#include "mgtaxa/exceptions.hpp"
#include "mgtaxa/debug.hpp"

#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>

namespace MGT {
	
inline int ipow(int base, int power) throw (std::domain_error) {
	if (power < 0) throw std::domain_error("power must be non-negative integer");
	if (power == 0) return 1;
	if (power == 1) return base;
	if (power % 2 == 0) return ipow(base * base, power / 2);
	if (power % 2 == 1) return base * ipow(base * base, power / 2);
}

/** @enum I_DEGEN Integer code for degenerate nucleotide.
It has to be zero. The enum mnemonic is created mainly for
better readability of the program source.*/

enum { I_DEGEN = 0 };	

/** Alphabet convertor from one-letter character to integer index.*/
 
class AbcConvCharToInt {

	public:
	
	enum { I_DEGEN = MGT::I_DEGEN };

	AbcConvCharToInt(const std::string& abc, const std::string& abcRevCompl);
	
	INuc operator()(CNuc c) const;
	INuc toINuc(CNuc c) const;
	CNuc toCNuc(INuc i) const;
	int nAbc() const;
	int nCodes() const;

	INuc revCompl(INuc i) const;	
	
	protected:

	/**If abc is 'ACTG', then m_CNucToINuc['A'] => 1, m_CNucToINuc['C'] => 2 and so on. 0 is reserved for 'N'*/
	CNuc m_CNucToINuc[g_maxCNuc];
	
	/**Extended alphabet, e.g. - 'NACGT'*/ 
	std::string m_abcExt;
	
	/**number of "real" alphabet symbols*/
	int m_nAbc;
	/**total number of encoded alphabet symbols (m_nAbc + 1)*/
	int m_nCodes;
	
	/**maps INuc index into reverse-complement INuc index.*/
	INuc m_iNucRevCompl[g_maxINuc];

};		


extern std::string g_defNucAbc;
extern std::string g_defNucAbcRevComp;
extern AbcConvCharToInt g_defAbcConvCharToInt;

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
	std::ostream& print(std::ostream& out, int kmerLen = g_maxKmerLen) const;
};

bool operator==(const Kmer& x, const Kmer& y);

std::ostream& operator<< (std::ostream& out, const Kmer& kmer);

std::string kmerToStr(const Kmer& kmer, int kmerLen, const AbcConvCharToInt& conv);

void kmerToRevCompl(const Kmer& kmer, int kmerLen, const AbcConvCharToInt& conv, Kmer& kmer);
  

/** Forward declaration of payload class for  KmerState.*/

class KmerStateData;

/** Pointer to KmerStateData.*/

typedef KmerStateData *PKmerStateData;

/** Type for unique ID of a k-mer. It should be stable between program invocations.
 * For dense  KmerStates implementation, this is defined as int.
 * For sparce (e.g. hash) representations, it might be k-mer itself.*/
 
typedef int KmerId;

/** Class that holds a counter for one k-mer along with jump pointers for next k-mers.*/

class KmerState {
	public:

	KmerState();

	typedef KmerState* PKmerState;
	
	std::ostream& print(std::ostream& out, const PKmerState pFirstState) const;

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

	const PKmerState getState() const {
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
 * Any implementation of KmerStates must guarantee this:
 * All degenerate states point to the very first (zero) state as their reverse-complement state.
 * That allows KmerCounter to link zero state to a single dummy KmerStateData object.
 * Thus, reverse-complement link is not symmetric for degenerate states.
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
		KmerStates(int kmerLen, const AbcConvCharToInt *pAbcConv);
		PKmerState nextState(PKmerState pState,INuc c);
		PKmerState revCompState(PKmerState pState);
		PKmerState revCompStateFirst(PKmerState pState);		
		bool isRevComp(PKmerState pState) const;
		bool isDegenState(PKmerState pState) const;
		PKmerState firstState();
		const PKmerState firstState() const;
		int numStates() const;
		KmerId idState(PKmerState pState) const;
		PKmerStateData getData(PKmerState pState);
		void setData(PKmerState pState,PKmerStateData pStateData);
		const Kmer& kmerState(PKmerState pState) const;
		std::ostream& print(std::ostream& out) const;
	private:
	/** Private copy constructor makes instances non-copyable.*/
    KmerStates( const KmerStates& );
    /** Private assignment operator makes instances non-copyable.*/
    const KmerStates& operator=( const KmerStates& );
	protected:
		PKmerState kmerToState(const Kmer& kmer) const; 
		void nextKmer(const Kmer& currKmer, INuc c, Kmer& nextKmer) const;
		int kmerToIndex(const Kmer& kmer) const;
		void indexToKmer(int ind,Kmer& kmer) const;
		void initAllKmers();
		void initRevCompl();		

	protected:
	const AbcConvCharToInt *m_pAbcConv;	
	int m_kmerLen;
	std::vector<int> m_blockSize;
	std::vector<int> m_blockStart;
	
	StateArray m_states;
	KmerArray m_kmers;
	
	/** Index of the first non-degenerate state in m_states.*/ 
	int iFirstStateNonDegen;
	
};


/** Class that computes k-mer counts of incoming nucleotide sequence.
It processes the input sequence and extracts the counts.
Internally, it maintains the KmerStateData payload data and moves through
states of the KmerStates state machine in response to incoming nucleotides.*/
 
class KmerCounter {
	
	public:
	
	KmerCounter(int kmerLen, const AbcConvCharToInt *pAbcConv = 0);
	~KmerCounter();
	
	void doCNuc(CNuc cnuc);
	void doINuc(INuc inuc);

	/** @name Interface to extract the results.
 	* The following set of methods defines the result extraction protocol.
 	* They must be called in a strict order because they change the internal state
 	* of KmerCounter object. The interface is designed for efficiency and flexibility
 	* of the caller.
 	* After doCNuc() has been called any number of times (we call this "accumulation cycle", 
 	* the extraction is done as in the
 	* following code sample, after which doCNuc() can be called again to accumulate new counts.
 	* @code
 	* int n = o.numKmers();
 	* o.startKmer();
 	* for(int i = 0; i < n; i++,o.nextKmer()) {
 	*     cout << o.getKmerId() << ":" << o.getKmerCount() << "\n";
 	* }
 	* o.finishKmer();
 	* @endcode
 	* nextKmer() can be called less than numKmers() times.
 	* */
	/*@{*/ 	
	/** Return number of k-mers found so far in current accumulation cycle.*/
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
	
	private:
	/** Private copy constructor makes instances non-copyable.*/
    KmerCounter( const KmerCounter& );
    /** Private assignment operator makes instances non-copyable.*/
    const KmerCounter& operator=( const KmerCounter& );
    
	protected:

	const AbcConvCharToInt *m_pAbcConv;
	
	PKmerState m_pSt;
	KmerStates *m_pStates;
	
	int m_kmerLen;
	
	/** Preallocated array of  KmerStateData objects.
	 * Array size is equal to the total number of states.
	 * This guarantees that there is always enough data elements
	 * regardless of the input sequence length.
	 * @todo Current array size is more than 50% excessive because
	 * we only store counts for one of every two reverse-complement
	 * states, and all degenerate states have zero state set as
	 * their reverse complement. This is a low priority optimization.*/
	std::vector<KmerStateData> m_data;
	
	/** One dummy KmerStateData object that is linked to the zero state,
	 * which in turn is set as a reverse-complement one for all other
	 * degenerate states. This way, it serves as a sink counter for
	 * all degenerate states. That in turn removes one branch condition
	 * from the time-critical code in doINuc().*/ 
	KmerStateData m_dataDegen;
	
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

/** Link given KmerState object to the given KmerStateData object.
 * This sets only uni-directed link from KmerState to KmerStateData.
 * @param pState - pointer to state
 * @param pStateData - pointer to data
 */   

inline void KmerStates::setData(PKmerState pState,PKmerStateData pStateData) {
	pState->m_pData = pStateData;
}

/** Return reference to Kmer object for a given KmerState obejct.
 * @param pState - pointer to state
 * This is done in constant time based on position of pState in internal array.
 */ 

inline const Kmer& KmerStates::kmerState(PKmerState pState) const {
	return m_kmers[pState-firstState()];
}

/** Return unique ID of the state's k-mer.
 * This is stable between invocations of the program.*/

inline KmerId KmerStates::idState(PKmerState pState) const {
	return pState->m_id;
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
 */ 

inline int KmerStates::kmerToIndex(const Kmer& kmer) const {
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

inline void KmerStates::indexToKmer(int ind, Kmer& kmer) const {
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
 * @todo make sure that reverse-complement pointers work correctly
 * for degenerate k-mer states and palyndromic states.
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
	return m_pStates->idState(m_data[m_iDataExtr].getState());
}
 
inline std::string KmerCounter::getKmerStr() const {
	return kmerToStr(m_pStates->kmerState(m_data[m_iDataExtr].getState()), m_kmerLen, *m_pAbcConv);
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


} // namespace MGT

#endif /*MGT_KMERS_HPP_*/
