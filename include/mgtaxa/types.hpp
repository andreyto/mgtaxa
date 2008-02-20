#ifndef MGT_TYPES_HPP_
#define MGT_TYPES_HPP_

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

} // namespace

#endif /*MGT_TYPES_HPP_*/
