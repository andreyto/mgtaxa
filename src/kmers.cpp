//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the MGTAXA package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#include "mgtaxa/kmers.hpp"

namespace MGT {

/** Print KmerState object for debugging.
@param pFirstState - pointer to the first KmerState, used here to convert internal pointers to indices
*/

std::ostream& KmerState::print(std::ostream& out, const PKmerState pFirstState) const {
	for(int i = 0; i < g_maxINuc; i++) {
		out << (m_next[i] - pFirstState) << '\t';
	}
	out << m_pData << '\t';
	out << (m_revComp - pFirstState) << '\t';
	out << m_isRevComp << '\t' << m_id << '\n';
	return out;
}

/** Constructor.
 * @param kmerLen - length of k-mers
 * @param pAbcConv - pointer to AbcConvCharToInt alphabet convertor object 
 * (stored inside this KmerStates object but not managed).
 * @param firstIdState - start k-mer ids from this value.
 */

KmerStates::KmerStates(int kmerLen, 
                       const AbcConvCharToInt  *pAbcConv,
                       RC_POLICY revCompPolicy,
                       KmerId firstIdState) {
    m_revCompPolicy = revCompPolicy;
    m_firstIdState = firstIdState;
    m_lastIdState = m_firstIdState; // will be computed during initialization
	m_pAbcConv = pAbcConv;
	int nAbc = m_pAbcConv->nAbc();
	if( kmerLen > g_maxKmerLen ) {
		throw KmerErrorLimits("kmerLen is too large");
	}
	if( nAbc + 1 > g_maxINuc ) {
		throw KmerErrorLimits("nAbc is too large");
	}
	m_kmerLen = kmerLen;
	m_blockSize.resize(kmerLen+1);
	m_blockStart.resize(kmerLen+1);
	int nAbcKmers = ipow(nAbc,kmerLen);
	int nKmers = nAbcKmers;
	for(int i = kmerLen, nPreKmers = nAbcKmers; i >= 0; i--) {
		m_blockSize[i] = nPreKmers;
		nPreKmers /= nAbc;		
		nKmers += nPreKmers;
	}
	m_blockStart[0] = 0;
	for(int i = 1; i <= kmerLen; i++) {
		m_blockStart[i] = m_blockSize[i-1] + m_blockStart[i-1];
	}
	m_states.resize(nKmers);
	m_kmers.resize(nKmers);
	initAllKmers();
	initRevCompl();
}


/** Initialize all KmerState and Kmer objects.*/

void KmerStates::initAllKmers() {
	int nCodes = m_pAbcConv->nCodes();	
	PKmerState pStateFirst = &m_states[0];
	for(int iState = 0; iState < m_states.size(); iState++) {
		PKmerState pState = &m_states[iState];
		Kmer& cKmer = m_kmers[iState];
#if ADT_DBG_LEVEL >= 4
		Kmer cKmerBak = cKmer;
#endif
		indexToKmer(iState,cKmer);
#if ADT_DBG_LEVEL >= 4
		int iStateCheck = kmerToIndex(cKmer);
		//cKmer.print(std::cerr) << "iState = " << iState << '\n';
		assert(iStateCheck == iState);
		//if the kMer was already non-zero, we assert that
		//the last call to indexToKmer did not change it:
		bool isNonZeroKmer = false;
		for(int i = 0; i < m_kmerLen; i++) {
			if( cKmerBak[i] != I_DEGEN ) {
				isNonZeroKmer = true;
				break;
			}
		}
		if( isNonZeroKmer ) {
			if(! (cKmer == cKmerBak) ) {
				ADT_LOG << "Previously filled k-mer content does not match: " <<
					ADT_OUTVAR(cKmer) << "  :  " << ADT_OUTVAR(cKmerBak) << "\n";
				ADT_ALWAYS(0);
			}
		}
#endif
		for(INuc c = 0; c < nCodes; c++) {
			Kmer nKmer;
            //**@todo Maybe use std::next_permutation()*/ 
			nextKmer(cKmer,c,nKmer);
			PKmerState pStateNext = kmerToState(nKmer);
			pState->m_next[c] = pStateNext;
			if( pStateNext > pState ) {
				// all omitted kmers (like a0bc and not like 00ab or abcd)
				// will have pStateNext pointing to the 0 position State,
				// so everything else must be stored to be processed in some next
				// iState loop iteration
				m_kmers[pStateNext - pStateFirst] = nKmer;
			}
		}
	}
}


/** Initialize reverse-complement links in KmerState objects.*/

void KmerStates::initRevCompl() {
	int idState = m_firstIdState;
	PKmerState pStateFirst = &m_states[0];
	for(int iState = 0; iState < m_states.size(); iState++) {
		PKmerState pState = &m_states[iState];
		Kmer& cKmer = m_kmers[iState];
		if( isDegenState(pState) ) {
			pState->m_revComp = pStateFirst;
			pState->m_isRevComp = true;
		}
		else {
			Kmer rKmer;
			kmerToRevCompl(cKmer, m_kmerLen, *m_pAbcConv, rKmer);
			PKmerState pStateRC = kmerToState(rKmer);
			// We intentionally do a redundant pass through all rev-comp
			// pairs to use the following asserts as a sanity check on our
			// k-mer calculation code.
			if( pState->m_revComp ){
				assert(pState->m_revComp == pStateRC);
			}
			if( pStateRC->m_revComp ){
				assert(pStateRC->m_revComp == pState);
			} 
			pState->m_revComp = pStateRC;
			pStateRC->m_revComp = pState;
            if( m_revCompPolicy == RC_MERGE )
            {
                if( pState <= pStateRC ) {
                    // In that order to work for palyndromes
                    pStateRC->m_isRevComp = true;
                    pState->m_isRevComp = false;
                }
                else {
                    pState->m_isRevComp = true;
                    pStateRC->m_isRevComp = false;
                }
            }
            else {
                pState->m_isRevComp = false;
            }
            if( ! pState->m_isRevComp ) {
				// All other m_id are left to be 0 because we
				// should never see them in the output
				pState->m_id = idState++;
            }
		}
	}
	pStateFirst->m_isRevComp = false;
    m_lastIdState = idState;
}

/** Print KmerStates object for debugging.*/

std::ostream& KmerStates::print(std::ostream& out) const {
	out << ADT_OUTVAR(m_kmerLen) << ADT_OUTVAR(m_states.size()) << "\n";
	PKmerState pFirstState = firstState();
	for(int iState = 0; iState < m_states.size(); iState++) {
		const PKmerState pState = const_cast<const PKmerState>(&m_states[iState]);
		const Kmer& cKmer = m_kmers[iState];
		out << ADT_OUTVAR(iState);
		out << "State: " << pState->print(out,pFirstState) << '\t' << ADT_OUTVAR(cKmer) 
				<< kmerToStr(cKmer,m_kmerLen,*m_pAbcConv) << "\n";
	}
	return out;
}

/**A constructor.
 * @param abc - a sequence of allowed non-degenerate character alphabet symbols (such as ACGT),
 * anything else that will be seen in the future sequence input
 * will be treated as degenerate symbols, equal to each other.
 * The index representation of degenerates is always 0.
 * @param abcRevCompl - a sequence of reverse-complement symbols for each element of abc (such as TGCA).
*/

AbcConvCharToInt::AbcConvCharToInt(const std::string& abc, const std::string& abcRevCompl) {
	m_abcExt = "N"+abc;
	int nAbc = abc.size();
	if( nAbc > g_maxINuc - 1 ) {
		throw KmerBadAlphabet();
	}
	if( abc.size() != abcRevCompl.size() ) {
		throw KmerBadAlphabet();
	}

	for(int i = 0; i < g_maxCNuc; i++) {
		m_CNucToINuc[i] = I_DEGEN;
	}
	for(int i = 0, ind = 1; i < nAbc; i++) {
		CNuc cnuc = abc[i];
		if( int(cnuc) > g_maxCNuc || cnuc < 1 ) {
			throw KmerBadNuc(cnuc);
		}
		m_CNucToINuc[cnuc] = ind++;
	}
	m_nAbc = nAbc;
	m_nCodes = nAbc + 1;
	
	// Init rev-compl index map
	
	for(int i = 0; i < g_maxINuc; i++) { 
		m_iNucRevCompl[i] = I_DEGEN;
	}

	for(int i = 0; i < nAbc; i++) {
		CNuc cnuc = abc[i];
		CNuc rcnuc = abcRevCompl[i];
		m_iNucRevCompl[toINuc(cnuc)] = toINuc(rcnuc);
	}
}


/** Default value for nucleotide alphabet */
std::string g_defNucAbc = "ACGT";

/** Default value for reverse complement of nucleotide alphabet */
std::string g_defNucAbcRevCompl = "TGCA";

AbcConvCharToInt g_defAbcConvCharToInt(g_defNucAbc,g_defNucAbcRevCompl);

/** A constructor.
 * @param kmerLen is a length of a k-mer. In the current implementation, all kmers are
 * precalculated and stored in memory, so be reasonable with this parameter.
 * @param pAbcConv is a to AbcConvCharToInt alphabet convertor
 * (stored inside this KmerCounter object but not managed).
 * @param revCompPolicy what to do about reverse complement k-mers
 * @param firstIdState Start k-mer IDs from this value (default 1 as in SVMLight)
*/

KmerCounter::KmerCounter(int kmerLen, 
                         const AbcConvCharToInt  *pAbcConv,
                         RC_POLICY revCompPolicy,
                         KmerId firstIdState) {
	if( pAbcConv == 0 ) {
		pAbcConv = & g_defAbcConvCharToInt;
	}
	m_pAbcConv = pAbcConv;
	m_kmerLen = kmerLen;
    m_revCompPolicy = revCompPolicy;
	m_pStates = new KmerStates(kmerLen,pAbcConv,m_revCompPolicy,firstIdState);
	m_data.resize(m_pStates->numStates());
	PKmerState pStateZero = m_pStates->firstState();
	m_dataDegen.setState(pStateZero);
	m_pStates->setData(pStateZero,&m_dataDegen);
	m_iDataEnd = 0;
	m_iDataExtr = 0;
	//set current state to be the first state
	m_pSt = m_pStates->firstState();
}

KmerCounter::~KmerCounter() {
	delete m_pStates;	
}


/** A constructor.
 * @param kmerLen is a max length of a k-mer. K-mers from 1 to k will be counted.
 * In the current implementation, all kmers are
 * precalculated and stored in memory, so be reasonable with this parameter.
 * @param pAbcConv is an alphabet convertor
 * (stored inside this KmerCounter object but not managed).
 * @param revCompPolicy what to do about reverse complement k-mers (currently only
 * RC_DIRECT is supported).
*/
KmerCounterLadder::KmerCounterLadder(int kmerLen, 
                         const AbcConvCharToInt  *pAbcConv,
                         RC_POLICY revCompPolicy) {
	if( pAbcConv == 0 ) {
		pAbcConv = & g_defAbcConvCharToInt;
	}
    ADT_ALWAYS(revCompPolicy == RC_DIRECT);
	m_pAbcConv = pAbcConv;
	m_kmerLen = kmerLen;
    m_revCompPolicy = revCompPolicy;
    m_counters.resize(kmerLen+1); //will index by k as in k-mer
    KmerId firstIdState = 1;
    // it leads to a little more simple code down the road if we output
    // the higher order k-mers first in the combined feature vector,
    // because we have to clean the accumulated state machine objects in that order.
    // So we set the ID ranges in that order.
    for(int i=kmerLen; i>=1; i--) {
        m_counters[i].reset(new KmerCounter(i,pAbcConv,m_revCompPolicy,firstIdState));
        firstIdState=m_counters[i]->getLastIdState();
    }
    makeTopDownKmerInd();
}

/** Create internal indices that connect each top level k-mer with (k-1)-mer prefix and 1-mer suffix.
 * Must be called from the constructor.*/

void KmerCounterLadder::makeTopDownKmerInd() {
    m_topBotDep.resize(m_kmerLen+1);
    m_topOneDep.resize(m_kmerLen+1);
    for(int k=m_kmerLen; k>1; k--) {
        KmerCounter& top = *m_counters[k];
        KmerCounter& bot = *m_counters[k-1];
        KmerCounter& one = *m_counters[1];
        KmerStates& topStates = top.getStates();
        KmerStates& botStates = bot.getStates();
        KmerStates& oneStates = one.getStates();
        indvec& topBotDep = m_topBotDep[k];
        indvec& topOneDep = m_topOneDep[k];
        topBotDep.resize(topStates.numStates());
        topOneDep.resize(topStates.numStates());
        for(Ind topInd=0; topInd < topStates.numStates(); topInd++) {
            const Kmer& topKmer = topStates.getKmer(topInd);
		    if(! topStates.isDegenState(topStates.getState(topInd)) ) {
                Kmer botKmer = topKmer.subKmer(0,k-1);
                Kmer oneKmer = topKmer.subKmer(k-1,k);
                Ind botInd = botStates.kmerToIndex(botKmer);
                Ind oneInd = oneStates.kmerToIndex(oneKmer);
                topBotDep[topInd] = botInd;
                topOneDep[topInd] = oneInd;
            }
        }
    }
}

} // namespace MGT

