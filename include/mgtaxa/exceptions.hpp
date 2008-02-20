#ifndef MGT_EXCEPTIONS_HPP_
#define MGT_EXCEPTIONS_HPP_

#include "mgtaxa/types.hpp"

#include <exception>

namespace MGT {

/** Base Kmer exception class.*/

class KmerError: public std::exception
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

/** Kmer exception to raise when supplied alphabet string exceeds implementation datatype limits.*/

class KmerBadAlphabet: public KmerError
{
public:
	KmerBadAlphabet():
	KmerError("Alphabet is too long")
	{}
};

/** Kmer exception to raise when some internal limit is exceeded by supplied input data.*/

class KmerErrorLimits: public KmerError
{
public:
	KmerErrorLimits(const char* msg): KmerError(msg)
	{}
};


/** Kmer exception to raise when unallowed one-letter nucleotide code is received in the input. */  

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

} // namespace MGT

#endif /*MGT_EXCEPTIONS_HPP_*/
