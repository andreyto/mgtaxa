//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the MGTAXA package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
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

/** Upper bound on k-mer length.
 * Selected to limit the required memory to ~1.7G on x86_64*/

const int g_maxKmerLen = 12;

} // namespace

#endif /*MGT_TYPES_HPP_*/
