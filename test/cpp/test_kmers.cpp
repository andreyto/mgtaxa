#include "mgtaxa/kmers.hpp"
#include "mgtaxa/version.hpp"
#include <boost/test/minimal.hpp>
#include <map>

namespace {

class RevComplementor {

	protected:

	enum { TableSize = 256 };
	
	public:
	
	RevComplementor(bool fastInit=false);

	void operator() (std::string& s);
		
	protected:

	char table[TableSize];
};


RevComplementor::RevComplementor(bool fastInit) {
	if( ! fastInit ) {
		for(int i = 0; i < TableSize; i++) {
			table[i] = i;
		}
	}
	table['a'] = 't';
	table['c'] = 'g';
	table['t'] = 'a';
	table['g'] = 'c';
	table['A'] = 'T';
	table['C'] = 'G';
	table['T'] = 'A';
	table['G'] = 'C';
}

void RevComplementor::operator() (std::string& s) {
	int  i = 0, j = s.size() - 1;
	for  ( ;  i < j;  i ++, j --) {
		char ch = table[s [i]];
		s [i] = table[s [j]];
		s [j] = ch;
	}
	if  (i == j) {
		s [i] = table[s [i]];
	}
}

RevComplementor g_revComp;

typedef std::map<std::string,int> KmerStrMap;

class KmerCounterNaive {
public:


KmerCounterNaive() {}

void fillFromSeq(const std::string& seq, int kmerLen) {
	m_counts.clear();
	for(int i = 0; i < seq.size() - kmerLen + 1; i++) {
		std::string kmer = seq.substr(i,kmerLen);
		//ADT_LOG << ADT_OUTVAR(kmer);
		if( kmer.find_first_not_of("ACGT") == std::string::npos ) {
			m_counts[kmer]++;
			g_revComp(kmer);
			//ADT_LOG << ADT_OUTVAR(kmer) << "\n";
			m_counts[kmer]++;
		}
	}
}

void fillFromCounter(const std::string& seq, MGT::KmerCounter& counter) {
	m_counts.clear();	
	for(int i = 0; i < seq.size(); i++) {
		counter.doCNuc(seq[i]);
	}
	int n = counter.numKmers();
	counter.startKmer();
	for(int i = 0; i < n; i++,counter.nextKmer()) {
		std::string kmerStr = counter.getKmerStr();
		int cnt = counter.getKmerCount();
    	m_counts[kmerStr] = cnt;
    	std::string kmerStrRC = kmerStr;
    	g_revComp(kmerStrRC);
    	if(! (m_counts.find(kmerStrRC) == m_counts.end() || kmerStrRC == kmerStr) ) {
    		throw "Rev-compl k-mer seen in KmerCounter output";
    	}
    	m_counts[kmerStrRC] += cnt; //increment so that palyndromic k-mers (like AT) get double counts
	}
	counter.finishKmer();	
}

const KmerStrMap& getCounts() const {
	return m_counts;
}

std::ostream& print(std::ostream& out) const {
	for( KmerStrMap::const_iterator p = m_counts.begin(); p != m_counts.end(); p++)	 {
		out << '(' << p->first << ':' << std::setw(3) << p->second << ") ";
	}
	return out;
}

KmerStrMap m_counts;

};

void compareMaps(const KmerStrMap& x, const KmerStrMap& y, KmerStrMap& absent, KmerStrMap& different) {
	for(KmerStrMap::const_iterator p_x = x.begin(); p_x != x.end(); p_x++) {
		KmerStrMap::const_iterator p_y = y.find(p_x->first);
		if( p_y != y.end() ) {
			if( p_y->second != p_x->second ) {
				different[p_x->first] = p_x->second - p_y->second;
			}
		}
		else {
			absent[p_x->first] = p_x->second;
		}
	}
}

extern std::string g_testingSeq;

bool test_KmerCounterCtor() {
	MGT::KmerCounter counter(2);
	return true;
}

bool test_KmerCounterNaive() {
	int kmerLen = 8;
	std::string testSeq = g_testingSeq.substr(0,g_testingSeq.size());
	for(int i = 0; i < testSeq.size() - 2; i += 37) { testSeq[i] = 'N'; testSeq[i+1] = 'N'; }
	ADT_LOG << "testSeq : \n" << testSeq << "\n";
	std::string testSeqRC = testSeq;
	g_revComp(testSeqRC);
	ADT_LOG << "testSeqRC : \n" << testSeqRC << "\n";	
	KmerCounterNaive counterNaive;
	counterNaive.fillFromSeq(testSeq,kmerLen);
	//ADT_LOG << "counterNaive : \n";
	//counterNaive.print(ADT_LOG) << "\n";
	MGT::KmerCounter counter(kmerLen);
	KmerCounterNaive counterMgt;	
	counterMgt.fillFromCounter(testSeq,counter);
	//ADT_LOG << "counterMGT : \n";
	//counterMgt.print(ADT_LOG) << "\n";
	const KmerStrMap& cntNai = counterNaive.getCounts();
	const KmerStrMap& cntMgt = counterMgt.getCounts();
	KmerStrMap cntAbsNaiMgt, cntAbsMgtNai, cntDiff;
	compareMaps(cntNai,cntMgt,cntAbsNaiMgt, cntDiff);
	cntDiff.clear();
	compareMaps(cntMgt,cntNai,cntAbsMgtNai, cntDiff);
	bool ret = true;
	if( cntAbsNaiMgt.size() ) {
		ADT_LOG << ADT_OUTVAR(cntAbsNaiMgt.size()) << '\n';
		ret = false;
	}
	if( cntAbsMgtNai.size() ) {
		ADT_LOG << ADT_OUTVAR(cntAbsMgtNai.size()) << '\n';
		ret = false;
	}
	if( cntDiff.size() ) {
		ADT_LOG << ADT_OUTVAR(cntDiff.size()) << '\n';
		ret = false;
	}
	ADT_LOG << ADT_OUTVAR(cntNai.size()) << '\n';
	ADT_LOG << ADT_OUTVAR(cntMgt.size()) << '\n';	
	return ret;
}

} // namespace

int add( int i, int j ) { return i+j; }

int test_main( int, char *[] )             // note the name!
{
	///BOOST_CHECK(test_KmerCounterCtor());
	BOOST_CHECK(test_KmerCounterNaive());
    // six ways to detect and report the same error:
/*    BOOST_CHECK( add( 2,2 ) == 4 );        // #1 continues on error
    BOOST_REQUIRE( add( 2,2 ) == 5 );      // #2 throws on error
    if( add( 2,2 ) != 4 )
      BOOST_ERROR( "Ouch..." );            // #3 continues on error
    if( add( 2,2 ) != 4 )
      BOOST_FAIL( "Ouch..." );             // #4 throws on error
    if( add( 2,2 ) != 4 ) throw "Oops..."; // #5 throws on error

    return add( 2, 2 ) == 4 ? 0 : 1;       // #6 returns error code*/
}

namespace {
std::string g_testingSeq =
"GATCTTACCCTCAAATATTATTAATCTACGACTATAAAGCTAACAATCCTTAAAAAGCTTGTCAATTTTTACTCTAAATA"
"ATAAAATCATAATTAACCAAAAACAAACTCTAAGAAGGAGGAGGAAAAATGAGTAAAATGCGTATTGGTGTATTAACAGG"
"TGGAGGCGATTGCCCAGGTCTAAACCCAGCTATCCGTGGTATTGTCATGAGAGCATTAGATTATGGAGACGAAGTTATAG"
"GTTTGAAGTATGGATGGGCTGGTCTTCTTAAGGCAGATACTATGCCTTTATCCTTAGAAATGGTAGAAGATATTCTTGAA"
"ATCGGCGGAACAATTCTTGGAAGTTCTAGAACAAACCCATTCAAAAAAGAAGAAGATGTTCAAAAATGTGTTGAGAACTT"
"CAAAAAGTTAAACTTAGATGCCTTAATCGCCATAGGTGGAGAAGACACTCTAGGAGTTGCATCAAAATTTAGCAAACTTG"
"GTCTTCCAATGATCGGAGTTCCAAAAACTATTGATAAAGATTTAGAAGAAACTGACTATACTCTTGGATTTGACACTGCT"
"GTTGAAGTAGTGGTAGATGCAATTAAAAGACTTAGAGATACTGCAAGATCTCATGCAAGAGTTATCGTAGTAGAAATAAT"
"GGGAAGACATGCAGGATGGTTAGCTCTTTATGGTGGGCTTGCAGGAGGAGCAGATTATATCTTAATCCCTGAAGTAGAAC"
"CTAATCTTGAGGATCTTTACAATCACATAAGAAAACTATACGCAAGAGGAAGAAATCACGCAGTTGTAGCCATCGCTGAG"
"GGAGTACAACTACCAGGATTTACTTATCAAAAAGGACAAGAAGGAATGGTAGATGCCTTTGGTCACATTCGCTTAGGTGG"
"TGTAGGTAATGTACTAGCCGAAGAGATACAGAAGAACTTGGGAATTGAAACCAGAGCCGTAATCTTAAGCCACCTACAAA"
"GGGGAGGAAGTCCATCAATAAGAGATAGAATCATGGGGCTTCTCCTTGGTAAGAAGGCTGTAGACTTAGTACATGAAGGA"
"AAATCTGGATTATTTGTTGCTGTAAAAGGAAACGAATTAGTACCAGTAGATATAACTTTAATTGAAGGGAAAACAAAGAA"
"TGTTGATCCTGCCTTTTATGAGAGCGTAAAAACTTTCTTTAATAAGTAAGAAAACTCTCAAAAAGAGGAAGGCGTTTGCC"
"AAGCTTTGATCGATGAATCGTTAGCATTAGATCCGCAGGATCCGTCAACATTATTGCTAGTAGGTATGGATGCTTTCTTC"
"ACAGCTAAGTACCAAAAAGCTATCAATTCATGGCAAACCATTCTCGACAGTAACCGCCCAGATGTAGATCGCGCTGCACT"
"GATGAATGCGATTGAAACAGCTAACTTACGTATTCAAGCCGAAACCTCTGGCATGCCTAACGATGATGTACACAAGCAGG"
"CAAAAGTCACAAGTAAATCAGTGAGTATTGCGATTTCAATCTCACCAGAGCTTGCGGCGAAAGCTGGTCCAGATGACACT"
"ATCTTTGTATTCGCCCGTGCAACTGAAGGCCCTAAAGTACCGTTAGCAGCGACCAAGGTCAGCGCAAAATCCTTACCTGT"
"AACGGTGACATTAGATGATAGCACTGGCATGGGGGGCGATGTTAAATTAAGCCAAACAGCAAGTGTTGAAGTCATTGCGG"
"TCTTATCTAAGCATGGCAATATCAAACCACAAGCTGGTGATATCCAAGGAAAGATTAGCAAAGTGAATGTGGGTGAAACC"
"GCAAACTTAGTACTGGATACTCAAGTACAATAAATCATCTACAATAGCGCTAAGTCTTGGCTTAGCGCAGGTATAAAAAA"
"ACCGGCTCATGCCGGTTTTTTTATACAGTAAAATTACTTAGCTTTAGACATAAACTCAATCGCAGCTTTGTAATCTTCAT"
"CGGTACAATCTGTACACATACCACCTGGTGGCATTGCATTTAATCCAGTTTTTACACTTTTAACTAGATTATCAACACCT"
"TTAGCTAAACGTGGTTCCCAGTCAGCAGTATTGTGGGATTTAGGAGCACCAGCAACACCCATGCTATGGCATACGGTACA"
"TGCCTTATTATAAATAGCTTCAGCATCTTGAGCTGAAACGTTGACTGACATAGTCAAAGCGGCGACTGCAGTCATGGCTA"
"ACAGTTTTTTCATGTTCAATGTTCCTAATGCAAATACGTACAAAGCTCTCGCTGCCCTTTGAGGCAGACTCATTATAGCG"
"GACTAAGATTACTACAAACTTCAATAAATCTCATCTTTTTGTGATTAGAATCTCAGCATAAATGATATGAGGCGACAAAC"
"GTGACAAATTCGTTAATATATTAGCAAGTATTAAAAAAGTTCGATTTGTACGAGAAATATTGCGGTTATTGGAACCATAT"
"CAAATTGCTGACATCTAGTCCTAAGCAATTGTAAACCCCGCAATAACGTTTGTGTTAGAATCACTTGTGTTTCCTTCCGG"
"CTTGAGATAGAGGGTTCAGTGACAAATATAATTTCAGTAGACACACTGTTATCAGCGAGCAAGCTGACATGCATCCGCGA"
"AGAACGCATCCTGTTTGATGAATTTAAGTTTTTGAGATTAATGCGGGCGATATCGTCCCAGATTGAAGGCCCTAATGGTG"
"CTGGCAAAACCAGTTTATTACGTATTTTGGCAGGCTTATCTCGACCCTATGCAGGGCAAACCTTTTATGTCAATGAAGAC"
"ATCAACCGATGTCGCGATGAATATAATGAAGATCTGCTGTATTTAGGCCATCTTGCAGGGGTGAAATCTGAGTTAACTGC"
"TGAAGAAAACCTTAACTTCAATTTAAGAATCAGTGGCTATGATGATTTTGATACATCTGCCATATTGGCAAAAGTGAATT"
"TATCTGGATTTGAAGAAGCCCTTGCAGGGCATTTATCCGCAGGTCAGCATCGCAGAACAGCATTGGCAAGACTCTGGCAC"
"AATGATTGTAAAATATGGATCTTGGATGAACCTTTTACTGCGATAGATAAAAGAGGTGTCGAAGAACTTGAACAATTGTT"
"TATTAAACATGCGGATAATGGTGGCTGCGTGATCCTTACCACTCACCAAGATATGGGCATTATCAAAGATGATAGGCTTC"
"GTAAAATTCGTCTAGATTATCGCTTCGTATAAGGTTTGGCATAATGAAAAGAGGCATCAGCTTTACTCAAGCGTTTTTTA"
"CATTACTGCAGCGGGATCTAAAAATCGCGATTCGCCACCGCGGTGATATTTTTAACCCATTGTTATTCTTTATTATGGTT"
"GTCACCCTATTTCCACTTGGTATTGGTCCAGAACCGCAAATGTTGGCGCGTATCGCCCCAGGGATTATTTGGGTTGCGGC"
"ATTGCTCGCGTCGATGCTATCGCTTGAGCGCCTTTTTAAAGCCGATTTTAGTGATGGCAGCTTAGAGCAGATGCTGCTTA"
"GTCCACAACCGCTCGCCATTTTAGTATTGGCAAAAGTATTAGCTCACTGGTTACTCACGGGTGTCCCACTCATCCTTATC";
} // namespace
