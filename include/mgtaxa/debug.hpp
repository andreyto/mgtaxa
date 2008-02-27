#ifndef ADT_DEBUG_HPP_
#define ADT_DEBUG_HPP_

#include <iostream>
#include <cassert>

#if !defined(ADT_DBG_LEVEL)
#  define ADT_DBG_LEVEL 8
#endif

#define ADT_LOG std::cerr

#define ADT_ALWAYS(x) assert(x)

#define ADT_OUTVAR(x) #x << ":\t" << (x) << '\t'

#endif /*ADT_DEBUG_HPP_*/
