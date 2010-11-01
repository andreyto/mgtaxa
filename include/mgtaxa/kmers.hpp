//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the MGTAXA package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
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

#include <boost/shared_ptr.hpp>

namespace MGT {

/** Array index type */
typedef int Ind;

int ipow(int base, int power) throw (std::domain_error);

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
    inline void set(int i,INuc v) { m_data[i] = v; }
    inline INuc get(int i) const { return m_data[i]; }
    inline Kmer subKmer(int first,int last) const { 
        Kmer x;
        for(int j=0; first < last; first++,j++) {
            x.set(j,get(first));
        }
        return x;
    }
	std::ostream& print(std::ostream& out, int kmerLen = g_maxKmerLen) const;
};

bool operator==(const Kmer& x, const Kmer& y);

std::ostream& operator<< (std::ostream& out, const Kmer& kmer);

std::string kmerToStr(const Kmer& kmer, int kmerLen, const AbcConvCharToInt& conv);

void kmerToRevCompl(const Kmer& kmer, int kmerLen, const AbcConvCharToInt& conv, Kmer& kmerOut);
  

/** Forward declaration of payload class for  KmerState.*/

class KmerStateData;

/** Pointer to KmerStateData.*/

typedef KmerStateData *PKmerStateData;

/** Type for unique ID of a k-mer. It should be stable between program invocations.
 * For dense  KmerStates implementation, this is defined as int.
 * For sparse (e.g. hash) representations, it might be k-mer itself.*/
 
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
	count(0),
	m_pState(0), 
    m_id(0)
	{}
	
	PKmerState getState() {
		return m_pState;
	}

	const PKmerState getState() const {
		return m_pState;
	}
	
	void setState(PKmerState pState) {
		m_pState = pState;
        ADT_ALWAYS(pState);
        m_id = pState->m_id;
	}

    KmerId idState() const {
        return m_id;
    }

	/** k-mer count.*/
	ULong count;
	
	protected:
	
	/** Pointer to the owning  KmerState object.*/
	PKmerState m_pState;

    /** ID of k-mer, cached here from m_pState to avoid pointer dereferencing during sorting.*/
    KmerId m_id;
};

/** How to treat reverse complements */
enum RC_POLICY { 
    RC_MERGE, /** Sum up counts into the lexicographically first k-mer in each rev-comp pair */
    RC_DIRECT /** Treat each k-mer independently */
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
 * length and N is the number of possible k-mers + O(min(L,N)*log(min(L,N))) if results are sorted
 * by k-mer ID for extraction.*/ 

class KmerStates {
	protected:
	typedef std::vector<KmerState> StateArray;
	typedef std::vector<Kmer> KmerArray; 
	public:
		KmerStates(int kmerLen, 
                const AbcConvCharToInt *pAbcConv,
                RC_POLICY revCompPolicy=RC_MERGE,
                KmerId firstIdState=1);
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
		PKmerStateData getData(Ind ind);
		void setData(PKmerState pState,PKmerStateData pStateData);
		const Kmer& kmerState(PKmerState pState) const;
        Ind indState(const PKmerState pState) const;
        const PKmerState getState(Ind ind) const;
        const Kmer& getKmer(Ind ind) const;
        KmerId getFirstIdState() const;
        KmerId getLastIdState() const;
        int getNumIds() const;
		std::ostream& print(std::ostream& out) const;
	private:
	/** Private copy constructor makes instances non-copyable.*/
    KmerStates( const KmerStates& );
    /** Private assignment operator makes instances non-copyable.*/
    const KmerStates& operator=( const KmerStates& );
    public:
	Ind kmerToIndex(const Kmer& kmer) const;
	protected:
		PKmerState kmerToState(const Kmer& kmer) const; 
		void nextKmer(const Kmer& currKmer, INuc c, Kmer& nextKmer) const;
		void indexToKmer(int Ind,Kmer& kmer) const;
		void initAllKmers();
		void initRevCompl();		

	protected:
	const AbcConvCharToInt *m_pAbcConv;	
	int m_kmerLen;
    RC_POLICY m_revCompPolicy;
	std::vector<int> m_blockSize;
	std::vector<int> m_blockStart;
	
	StateArray m_states;
	KmerArray m_kmers;
	
	/** Index of the first non-degenerate state in m_states.*/ 
	int iFirstStateNonDegen;

    /** IDs of output (non-degenerate, non-reverse-complement) k-mers fill 
     * a half open interval [m_firstIdState,m_lastIdState).*/
    KmerId m_firstIdState, m_lastIdState;
};


/** Class that computes k-mer counts of incoming nucleotide sequence.
It processes the input sequence and extracts the counts.
Internally, it maintains the KmerStateData payload data and moves through
states of the KmerStates state machine in response to incoming nucleotides.*/
 
class KmerCounter {
	public:
	
	KmerCounter(int kmerLen, 
            const AbcConvCharToInt *pAbcConv = 0,
            RC_POLICY revCompPolicy=RC_MERGE,
            KmerId firstIdState=1);
	~KmerCounter();
	
	void doCNuc(CNuc cnuc);
	void doINuc(INuc inuc);

	int maxNumKmers(ULong seqLen) const;

	/** @name Interface to extract the results.
 	* The following set of methods defines the result extraction protocol.
 	* They must be called in a strict order because they change the internal state
 	* of KmerCounter object. The interface is designed for efficiency and flexibility
 	* of the caller.
 	* After doCNuc() has been called any number of times (we call this "accumulation cycle"), 
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
    * sumKmerCounts() can be called outside of startKmer()...finishKmer() block
 	* */
	/*@{*/ 	
	/** Return number of k-mers found so far in current accumulation cycle.*/
	int numKmers() const;
	/** Return sum of non-degenerate k-mer counts found so far in current accumulation cycle.
     * Can be only called outside of startKmer()...finishKmer() block.
     * Complexity: linear in the number of found km-mers.*/
	ULong sumKmerCounts() const;
	/** Return sum of degenerate k-mer counts found so far in current accumulation cycle.
     * Can be only called outside of startKmer()...finishKmer() block.
     * Complexity: constant time*/
	ULong sumDegenKmerCounts() const;
	/** Prepare internal state for result extraction.
     * @param doSort - if true, the results will be sorted (complexity will be N*log(N) where
     * N is numKmers(). SVM sparse feature vector representation needs sorted results.*/
	void startKmer(bool doSort=true);
	/** Advance internal state to extract next k-mer results.*/
	void nextKmer();
	/** Accessor to get count value from the currently extracted k-mer.*/
	ULong getKmerCount() const;
	/** Accessor to get Id from the currently extracted k-mer.*/
	int getKmerId() const;
	/** Accessor to get k-mer string (such as 'ACCCT') for the currently extracted k-mer.*/ 
	std::string getKmerStr() const;
	/** Accessor to get index of the state for the currently extracted k-mer.
     * Should be used only by implementation-aware code such as KmerCounterLadder*/ 
    Ind indState() const;
	/** Accessor to get pointer to the state for the currently extracted k-mer.
     * Should be used only by implementation-aware code such as KmerCounterLadder*/ 
    const PKmerState getState() const;
	/** Finalize result extraction cycle - new series of doCNuc() calls can be done afterwards.*/
	void finishKmer();
	/*@}*/ 
	
    KmerId getFirstIdState() const;
    KmerId getLastIdState() const;
    int getNumIds() const;

    public:
    
    /** Get reference to internal KmerStates array.
     * Declared public only to be used by implementation-aware code such as KmerCounterLadder.*/
    inline KmerStates& getStates() {
        return *m_pStates;
    }

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

    /** How reverse complements are treated - one of RC_XXX.*/
    RC_POLICY m_revCompPolicy;
	
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
	 * from the time-critical code in doINuc().
     * We also query from it the total count of degenerate states in the
     * last accumulation cycle.*/ 
	KmerStateData m_dataDegen;
	
	/** Index of the first unused element in ::m_data.*/
	int m_iDataEnd;
	
	/** Index of m_data element that is currently being extracted.*/
	int m_iDataExtr; 
};


/** Class that computes k-mer frequences for k = 1...n at the same time.
 * It normalizes observed frequencies for a given k with expected frequencies
 * calculated based on observed frequencies for (k-1) and 1.
 * The implementation:
 * 1. Compute counts for k = 1...n. Currently we just do it independently for
 * each k through calls to doCNuc(). Future implementation may compute n-mers and then find (n-1)-mer counts
 * by incrementing by one the count for (n-1) prefix of each n-mer plus special handling of the right boundary 
 * n-mer (and repeating recursively for (n-2) and so on). That will be more efficient for long sequences.
 * 2. Compute expected counts for each n-mer found in (1).
 * 3. Store normalized frequences into the output array.*/

class KmerCounterLadder {
    public:

	KmerCounterLadder(int kmerLen, 
            const AbcConvCharToInt *pAbcConv = 0,
            RC_POLICY revCompPolicy=RC_DIRECT);

	void doCNuc(CNuc cnuc);

    /** Accumulate counts for a CNuc sequence represented by an interator range*/
    template<class IterCNuc>
    inline
    void process(IterCNuc first,IterCNuc last) {
        for( ; first != last; ++first ) {
            doCNuc(*first);
        }
    }

    template<class IterInd,class IterVal>
    void counts(IterVal valObs,IterVal valExp,IterInd ind,IterInd sizes);
    
    template<class IterInd>
    int numKmers(IterInd sizes) const;

    template<class IterInd>
    int maxNumKmers(ULong seqLen,IterInd sizes) const;
    
    int numSubFeatures() const {
        return m_kmerLen;
    }

    protected:

    void makeTopDownKmerInd();

    protected:

    typedef boost::shared_ptr<KmerCounter> PKmerCounter;
    typedef std::vector<PKmerCounter> PKmerCounterArray;

    typedef std::vector<Ind> indvec;

    /** KmerCounter pointer array [0,m_kmerLen], with elements
     * [1,m_kmerLen] initialized */
    PKmerCounterArray m_counters;
    /** m_topBotDep[k_mer_index]->(k-1)_mer_index
     * array [0,m_kmerLen], with elements
     * [1,m_kmerLen] initialized*/
    std::vector<indvec> m_topBotDep;
    /** m_topOneDep[k_mer_index]->(1)_mer_index
     * array [0,m_kmerLen], with elements
     * [1,m_kmerLen] initialized*/
    std::vector<indvec> m_topOneDep;

	const AbcConvCharToInt *m_pAbcConv;
	
	int m_kmerLen;

    /** How reverse complements are treated - one of RC_XXX.*/
    RC_POLICY m_revCompPolicy;
};

} // namespace MGT

#include "mgtaxa/kmers_inc.hpp"

#endif /*MGT_KMERS_HPP_*/
