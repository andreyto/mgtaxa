
KmerStates::KmerStates(int kmerLen, int nAbc) {
	if( kmerLen > maxKmerLen ) {
		throw KmerErrorLimits("kmerLen is too large");
	}
	if( nAbc + 1 > maxINuc ) {
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
			PKmerState pStateNext = kmerToPointer(nKmer);
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



// alphabet must list all "real" nucleotide symbols (e.g. ATGC)
// anything else that will be seen in the future sequence input
// will be treated as N, equal to each other (e.g. X will be considered an N too)

KmerCounter::KmerCounter(int kmerLen, const CNuc abc[], int nAbc) {
	ctorCNucToINuc(abc,nAbc);
	ctorKmerArray(kmerLen);
}

void KmerCounter::ctorCNucToINuc(const CNuc abc[], int nAbc) {
	if( nAbc > maxINuc - 1 ) {
		throw KmerBadAlphabet();
	}
	m_CNucToINuc = new CNuc[maxCNuc];
	for(int i = 0; i < maxCNuc; i++) {
		m_CNucToINuc[i] = 0;
	}
	for(int i = 0, ind = 1; i < nAbc; i++) {
		CNuc cnuc = abc[i];
		if( cnuc > maxCNuc || cnuc < 1 ) {
			throw KmerBadNuc(cnuc);
		}
		m_CNucToINuc[cnuc] = ind++;
	}
	m_nAbc = nAbc;
	m_nINuc = nAbc + 1;
}

void KmerCounter::ctorKmerArray(int kmerLen) {
	m_kmerLen = kmerLen;
	//n-dim array iterator ndim = kmerLen, extent = m_nAbc
	//easier to use inverted iterator: i goes along entire flat array,
	//and we calculate n-dim index for each point
	for(int i = arr.beginPos(); i < arr.endPos(); i++) {}
}
		

~KmerCounter::KmerCounter() {
	delete [] m_CNucToINuc;
}


} // namespace Phy

int main() {
vector<char> a(1000*1000*100);
a[0] = 3;
a[1] = 200;
for(int i = 2; i < a.size(); i++) {
	a[i] = (a[i-1] + a[i-2])/2;
}
if( a[a.size()-1] < 10 ) {
	cout << a[a[a.size()-1]] << "\n";
}
else {
	cout << a[a[a.size()-1]/2] << "\n";
}
return 0;
}
