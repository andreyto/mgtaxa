#include <vector>
#include <iostream>
#include <exception>
#include <cmath>

namespace Phy {

inline int ipow(int b, int p) {
	return int(std::pow(b,p));
}	

// small integer type to encode a nucledotide,
// e.g. 0,1,2,3,4 (0 - for N)

typedef char INuc;

const int maxINuc = 5;

// type for letter encoded nucleotide

typedef char CNuc;

const int maxCNuc = 256;

typedef unsigned long ULong;

const int maxKmerLen = 10;

class KmerError: public exception
{
public:
	KmerError():
	m_msg("KmerError")
	{}
	KmerError(const char* msg):
	m_msg(msg)
	{}

  virtual const char* what() const throw()
  {
    return m_msg;
  }
 protected:
 const char* m_msg;
};

class KmerBadAlphabet: public KmerError
{
public:
	KmerBadAlphabet():
	m_msg("Alphabet is too long")
	{}
};

class KmerErrorLimits: public KmerError
{
public:
	KmerErrorLimits(const char* msg): KmerError(msg)
	{}
};

class KmerBadNuc: public KmerError
{
  
  KmerBadNuc(CNuc cnuc) 
  {
  	m_msg[0] = cnuc;
  	m_msg[1] = '\0';
  }
  
  virtual const char* what() const throw()
  {
    return m_msg;
  }
  
  protected:
  
  char m_msg[2];
  
};	


class KmerCounter {
	public:
	KmerCounter():
	count(0) {
		for(int i = 0; i < maxINuc; i++) {
			m_next[i] = 0;
		}
	}
	typedef KmerCounter* PKmerCounter;
	PKmerCounter m_next[maxINuc];
	ULong count;	
};

typedef KmerCounter::PKmerCounter PKmerCounter;

class Kmer {
	public:
	Kmer() {
		for(int i = 0; i < maxKmerLen; i++) { m_data[i] = 0; }
	}
	INuc m_data[maxKmerLen];
	inline INuc& operator[](int i) { return m_data[i]; }
	inline const INuc& operator[](int i) const { return m_data[i]; }	
};

class KmerWalker {
	protected:
	typedef std::vector<KmerCounter> CounterArray;
	typedef std::vector<Kmer> KmerArray; 
	public:
	//1. transition rule: nextItem(currItem,cNuc)->nextItem
	//2. itemToPointer(item) -> pointer
	//3. initAllItems()
	public:
		KmerArray(int kmerLen,int nAbc);
		PKmerCounter kmerToPointer(const Kmer& kmer); 
		inline void nextKmer(const Kmer& currKmer, INuc c, Kmer& nextKmer);
	protected:
		int kmerToIndex(const Kmer& kmer);
		void indexToKmer(int ind,Kmer& kmer);
		void initAllKmers();

	protected:
	//number of "real" alphabet symbols
	int m_nAbc;
	//total number of encoded alphabet symbols
	//(m_nAbc + 1)
	int m_nCodes;
	int m_kmerLen;
	std::vector<int> m_blockSize;
	std::vector<int> m_blockStart;
	
	CounterArray m_counts;
	KmerArray m_kmers;
	
};


inline PKmerCounter KmerWalker::kmerToPointer(const Kmer& kmer) {
	int ind = kmerToIndex(kmer);
	return & m_counts[ind]; 
}

inline int KmerWalker::kmerToIndex(const Kmer& kmer) {
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


inline int KmerWalker::indexToKmer(int ind, Kmer& kmer) {
	int iBlock = 0;
	for(; iBlock < kmerLen; iBlock++) {
		if( m_blockStart[iBlock] > ind ) {
			break;
		}
	}
	iBlock -= 1;
	int ndim = iBlock; // iBlock == number of non-0 dimensions
	int indBlock = ind - m_blockStart[ndim];
	extent = m_nAbc;
	int stride = m_nAbc; 
	for( int i = m_kmerLen - ndim; i < m_kmerLen; i++ ) {
		kmer[i] = indBlock % stride;
		stride *= extent;
	}
}

KmerWalker::KmerWalker(int kmerLen, int nAbc) {
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
	m_counts.resize(nKmers);
	m_kmers.resize(nKmers);
}

void KmerWalker::initAllKmers() {
	PKmerCounter pCounterFirst = &m_counts[0];
	PKmerCounter pCounter = pCounterFirst;
	for(int iCounter = 0; iCounter < m_counts.size(); iCounter++) {
		PKmerCounter pCounter = &m_counts[iCounter];
		Kmer& cKmer = m_kmers[iCounter];
		indexToKmer(iCounter,cKmer);
		//TODO: make assertion here that indexToKmer returns the same as
		//the content of cKmer if it was initialized before
		for(INuc c = 0; c < m_nAbc; c++) {
			Kmer nKmer;
			nextKmer(cKmer,c,nKmer);
			PKmerCounter pCounterNext = kmerToPointer(nKmer);
			pCounter->m_next[c] = pCounterNext;
			if( pCounterNext > pCounter ) {
				// all kmers that are like a0bc and not like 00ab or abcd
				// will have pCounterNext pointing to the 0 position counter,
				// so everything else must be stored to be processed in some next
				// iCounter loop iteration
				m_kmers[pCounterNext - pCounterFirst] = nKmer;
			}
		}
	}
}

class KmerCounter {
	
	public:
	
	KmerCounter(int kmerLen);
	
	inline void nextCNuc(CNuc cnuc) {
		nextINuc(m_CNucToINuc[cnuc]);
	}
	
	inline void nextINuc(INuc inuc) {
		m_pKmer = m_pKmer->next[inuc];
		m_pKmer->count++;
	}
	
	protected:
	
	void ctorCNucToINuc(const CNuc Abc[], int nAbc);
	
	protected:
	
	KmerArray m_kmers;
	KmerArrayPointer m_pKmer;
	const CNuc * m_CNucToINuc;
	//number of "real" alphabet symbols
	int m_nAbc;
	//total number of encoded alphabet symbols
	//(m_nAbc + 1)
	int m_nINuc;
	int m_kmerLen;
	
};

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
