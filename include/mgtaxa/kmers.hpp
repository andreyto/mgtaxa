#ifndef MGT_KMERS_HPP_
#define MGT_KMERS_HPP_

#include "mgtaxa/exceptions.hpp"

#include <vector>
#include <iostream>
#include <cmath>

namespace MGT {


/** Small integer type to encode a nucledotide.
 * Example: 0,1,2,3,4 encodes NACTG (N - degenerate symbol.*/

typedef char INuc;

/** Upper bound on valid INuc values.*/

const int g_maxINuc = 5;

/** Type for one-letter encoded nucleotide.
 * Example 'NACTG'.*/

typedef char CNuc;

/** Upper bound on valid CNuc values.*/

const int g_maxCNuc = 256;

typedef unsigned long ULong;

/** Upper bound on k-mer length.*/

const int g_maxKmerLen = 10;


/** Class Alphabet convertor from one letter character to integer index.
 * @see AbcConvCharToInt() constructor.*/
 
class AbcConvCharToInt {

	public:

	AbcConvCharToInt(const CNuc abc[], int nAbc);
	~AbcConvCharToInt();
	
	INuc operator()(CNuc c);
	
	protected:

	/**If abc is 'ACTG', then m_CNucToINuc['A'] => 1, m_CNucToINuc['C'] => 2 and so on.*/
	const CNuc * m_CNucToINuc;
	
	/**number of "real" alphabet symbols*/
	int m_nAbc;
	/**total number of encoded alphabet symbols (m_nAbc + 1)*/
	int m_nINuc;

};		


/** Class that holds a counter for one k-mer along with jump pointers for next k-mers.*/

class KmerState {
	public:

	KmerState();

	typedef KmerState* PKmerState;

	public:

	/** Declared public only for access by KmerStates class.*/
	PKmerState m_next[g_maxINuc];

	/** Declared public only for access by KmerCounter class.*/
	ULong m_count;	
};

/** Pointer type for @see KmerCounter */

typedef KmerState::PKmerState PKmerState;


/** Class that describes k-mer in integer index representation.*/ 

class Kmer {
	public:
	Kmer() {
		for(int i = 0; i < g_maxKmerLen; i++) { m_data[i] = 0; }
	}
	INuc m_data[g_maxKmerLen];
	inline INuc& operator[](int i) { return m_data[i]; }
	inline const INuc& operator[](int i) const { return m_data[i]; }	
};

/** State machine for k-mers.
 * It is optimized for the following use case:
 * 1. relatively short k-mers (k < 10)
 * 2. k-mer counts will be computed for many short chunks of sequence
 * (reset counts, stream through sequence chunk, extract the counts, repeat).
 * Typical length of sequence chunk (~1000 nuc) is less than the number of
 * possible k-mers (~16,000 for k = 7).
 * 3. k-mer counts will be extracted only for non-redundant subset of 
 * reverse complement k-mers (~8,000 for k = 7).
 * @see nextState implements the key operation of moving to the next k-mer
 * from currently seen one.
 * In view of (2), */ 

class KmerStates {
	protected:
	typedef std::vector<KmerState> StateArray;
	typedef std::vector<Kmer> KmerArray; 
	public:
	//1. transition rule: nextItem(currItem,cNuc)->nextItem
	//2. itemToPointer(item) -> pointer
	//3. initAllItems()
	public:
		KmerStates(int kmerLen,int nAbc);
		PKmerState nextState(PKmerState pState,INuc c);
		PKmerState firstState();
	protected:
		PKmerState kmerState(const Kmer& kmer); 
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
 * In the future, it can be easily parameterized on different @see KmerStates
 * implementations.*/
 
class KmerCounter {
	
	public:
	
	KmerCounter(int kmerLen);
	
	inline void doCNuc(CNuc cnuc) {
		doINuc((*m_pAbcConv)(cnuc));
	}
	
	inline void doINuc(INuc inuc) {
		m_pSt = m_pStates->nextState(m_pSt,inuc);
		m_pSt->m_count++;
	}
	
	protected:

	AbcConvCharToInt * m_pAbcConv;
	
	PKmerState m_pSt;
	KmerStates *m_pStates;
	
	int m_kmerLen;
	
};


////////////////////////////////////////////////////////////////
///
/// Inline implementations
///
////////////////////////////////////////////////////////////////

inline INuc AbcConvCharToInt::operator()(CNuc c) {
	return m_CNucToINuc[c];
}


inline KmerState::KmerState():
	m_count(0) {
		for(int i = 0; i < g_maxINuc; i++) {
			m_next[i] = 0;
		}
	}

inline PKmerState KmerStates::nextState(PKmerState pState,INuc c) {
	return pState->m_next[c];
}

inline PKmerState KmerStates::firstState() {
	return &m_states[0];
}


inline PKmerState KmerStates::kmerState(const Kmer& kmer) {
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

} // namespace MGT

#endif /*MGT_KMERS_HPP_*/
