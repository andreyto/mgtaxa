#include "mgtaxa/kmers.hpp"

namespace MGT {

KmerStates::KmerStates(int kmerLen, int nAbc) {
	if( kmerLen > g_maxKmerLen ) {
		throw KmerErrorLimits("kmerLen is too large");
	}
	if( nAbc + 1 > g_maxINuc ) {
		throw KmerErrorLimits("nAbc is too large");
	}
	m_kmerLen = kmerLen;
	m_nAbc = nAbc;
	m_nCodes = nAbc + 1;
	m_blockSize.resize(kmerLen+1);
	m_blockStart.resize(kmerLen+1);
	int nAbcKmers = ipow(m_nAbc,kmerLen);
	int nKmers = nAbcKmers;
	for(int i = kmerLen, nPreKmers = nAbcKmers; i > 0; i--) {
		m_blockSize[i] = nPreKmers;
		nPreKmers /= m_nAbc;		
		nKmers += nPreKmers;
	}
	m_blockStart[0] = 0;
	for(int i = 1; i < kmerLen; i++) {
		m_blockStart[i] = m_blockSize[i-1] + m_blockStart[i-1];
	}
	m_states.resize(nKmers);
	m_kmers.resize(nKmers);
}

void KmerStates::initAllKmers() {
	PKmerState pStateFirst = &m_states[0];
	PKmerState pState = pStateFirst;
	for(int iState = 0; iState < m_states.size(); iState++) {
		PKmerState pState = &m_states[iState];
		Kmer& cKmer = m_kmers[iState];
		indexToKmer(iState,cKmer);
		//TODO: make assertion here that indexToKmer returns the same as
		//the content of cKmer if it was initialized before
		for(INuc c = 0; c < m_nAbc; c++) {
			Kmer nKmer;
			nextKmer(cKmer,c,nKmer);
			PKmerState pStateNext = kmerToState(nKmer);
			pState->m_next[c] = pStateNext;
			if( pStateNext > pState ) {
				// all kmers that are like a0bc and not like 00ab or abcd
				// will have pStateNext pointing to the 0 position State,
				// so everything else must be stored to be processed in some next
				// iState loop iteration
				m_kmers[pStateNext - pStateFirst] = nKmer;
			}
		}
	}
}


/**A constructor.
 * @param abc is a sequence of allowed non-degenerate character alphabet symbolsal (e.g. ATGC),
 * anything else that will be seen in the future sequence input
 * will be treated as degenerate symbols, equal to each other. 
*/

AbcConvCharToInt::AbcConvCharToInt(const std::string& abc) {
	m_abcExt = "N"+abc;
	int nAbc = abc.size();
	if( nAbc > g_maxINuc - 1 ) {
		throw KmerBadAlphabet();
	}
	m_CNucToINuc = new CNuc[g_maxCNuc];
	for(int i = 0; i < g_maxCNuc; i++) {
		m_CNucToINuc[i] = 0;
	}
	for(int i = 0, ind = 1; i < nAbc; i++) {
		CNuc cnuc = abc[i];
		if( cnuc > g_maxCNuc || cnuc < 1 ) {
			throw KmerBadNuc(cnuc);
		}
		m_CNucToINuc[cnuc] = ind++;
	}
	m_nAbc = nAbc;
	m_nINuc = nAbc + 1;
}

AbcConvCharToInt::~AbcConvCharToInt() {
	delete [] m_CNucToINuc;
}


/** A constructor.
 * @param kmerLen is a length of a k-mer. In the current implementation, all kmers are
 * precalculated and stored in memory, so be reasonable with this parameter.
 * @param abc is a parameter for AbcConvCharToInt::AbcConvCharToInt() constructor.
*/

KmerCounter::KmerCounter(int kmerLen, const std::string& abc) {
	m_pAbcConv = new AbcConvCharToInt(abc);
	ctorKmerArray(kmerLen);
	m_data.resize(m_pStates->numStates())
	m_iDataEnd = 0;
	m_iDataExtr = 0;
}

void KmerCounter::ctorKmerArray(int kmerLen) {
	m_kmerLen = kmerLen;
	//n-dim array iterator ndim = kmerLen, extent = m_nAbc
	//easier to use inverted iterator: i goes along entire flat array,
	//and we calculate n-dim index for each point
	for(int i = arr.beginPos(); i < arr.endPos(); i++) {}
}
		

KmerCounter::~KmerCounter() {
	delete m_abcConv;
}


} // namespace MGT
