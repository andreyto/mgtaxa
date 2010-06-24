//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the MGTAXA package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//

#include <algorithm>
#include <iterator>
  /**
   *  @brief Remove consecutive values from a sequence using a predicate and store counts.
   *  @param  first        A forward iterator.
   *  @param  last         A forward iterator.
   *  @param  binary_pred  A binary predicate.
   *  @param  result       An output iterator.
   *  @return  An iterator designating the end of the resulting sequence in input range.
   *
   *  Removes all but the first element from each group of consecutive
   *  values for which @p binary_pred returns true.
   *  Counts for each group are stored in the output iterator range.
   *  unique() is stable, so the relative order of elements that are
   *  not removed is unchanged.
   *  Elements between the end of the resulting sequence and @p last
   *  are still present, but their value is unspecified.
   *  This is a modified STL algorithm from gcc 4.10
  */
  template<typename _ForwardIterator, typename _BinaryPredicate, typename _OutputIterator>
    _ForwardIterator
    unique_counts(_ForwardIterator __first, _ForwardIterator __last,
           _BinaryPredicate __binary_pred, _OutputIterator __result)
    {
        if (__first == __last)
	        return __last;
        typename std::iterator_traits<_ForwardIterator>::difference_type c=1;
        _ForwardIterator __dest = __first;
        while (++__first != __last)
        {
        	if (!__binary_pred(*__dest, *__first))
            {
            	*++__dest = *__first;
                *__result = c;
                ++__result;
                c = 1;
            }
            else
            {
                c++;
            }
        }
        *__result = c;
        ++__result;
        return ++__dest;
    }


#include <iostream>
#include <vector>
#include <string>

using namespace std;

int main() {
    string a = "aaaaaaabcdeeeefffgijkkkllll";
    //string a = "abcdeeeefffgijkkkl";
    //string a = "aaaa";
    cout << a << "\n";
    vector<int> vc(a.size());
    string::iterator last = unique_counts(a.begin(),a.end(),equal_to<char>(),vc.begin());
    copy(vc.begin(),last-a.begin()+vc.begin(),ostream_iterator<int>(cout, ""));
    cout << "\n";
    copy(a.begin(),last,ostream_iterator<char>(cout, ""));
    cout << "\n";
    return 0;
}

