//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the MGTAXA package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef ARRAY_BOOST_PY_H_
# define ARRAY_BOOST_PY_H_

/*
   Copyright: Andrei Tovchigrechko, 2001 

   The wrapper class around NumPy Numeric array object for boost::python library.

   The goal is to provide the object that is easy to pass to/from boost::python
   extensions and easy to extract the pointer to C data array and data layout. There are several good
   multidimensional array and matrix packages (Blitz, MTL, Pooma), which use same data storage
   scheme as NumPy (flat data array combined with 'dimensions' and 'strides' arrays). They allow
   for easy construction of their array objects around externally supplied data.
   Blitz and Pooma  use expression templates to avoid creation of temporaries in
   array ariphmetics and provide other compile-time optimizations. They also require type of data
   elements and rank to be specified as template parameters. It seems wise to use one of these 
   packages for actual array processing, rather than try to implement full NumPy interface in C++.
   The downside is that the initial data must be allocated as NumPy array object if it has to be
   returned into Python script without making the full copy, since Blitz and others do their own
   memory management with arrays that they allocate.

   The code also provides a set of conversion functions from other Python objects to 'array' objects 
   implemented as thin wrappers around NumPy functions. The functions are named 'to_array_xxx()', where
   suffix 'xxx' is a combination of (copy | contiguous | downcast). The suffix determines what function does if 
   argument is already a NumPy array, i.e. always return an array with non-shared copy of data or continuously allocated 
   data, allow for casting of data with loss of accuracy (e.g. returning array of 'float' if argument array has data
   in 'double').

   When extracting data pointer or making conversion to array object, run-time checks are performed 
   to ensure compatibility of data type stored in NumPy array with data type expected by C++ code.

   Problems: no matter how I import and imitialize Numeric module from within my C++ extension module (see
   example in 'test_array_ext.cpp'), I still have to do 'import Numeric' from Python before using my
   extension module or face the error message about import failure (Python 2.1.1).
   
*/

#include <complex>

# include "boost/python/objects.hpp"

// this is to avoid multiple definition of PyArray_API global
// variable, see "Numeric/arrayobject.h"
#define NO_IMPORT_ARRAY

#include "Numeric/arrayobject.h"

namespace boost { namespace python {


  // namespace to put stuff used for internal needs of array class

  namespace pyarray {

    // class to return run-time info (e.g. size of data element) for each PyArray_XXX
    // type constant
    
    struct type_info_rec {
      int itype; // value of PyArray_XXX constant 
      std::size_t size; //size of corresponding c data type
      const char *s_type; //itype name as string
    };
    
    class type_info {

      static const type_info_rec atype_info[];

    public:
      
      // return size of element with type code itype
      static std::size_t get_size(int itype);

      // return string name of type code itype
      static const char* get_s_type(int itype);
    };


    // the following classes establish compile-time
    // correspondence between C types and PyArray element
    // types, declared as enum constants in "Numeric/arrayobject.h"

    template<class c_type> struct c_to_py {};

    template<> struct c_to_py<char> 
    { enum { py_type = PyArray_CHAR }; };

    template<> struct c_to_py<unsigned char> 
    { enum { py_type = PyArray_UBYTE }; };

    // we omit SBYTE since 'signed char'=='char'
    template<> struct c_to_py<short> 
    { enum { py_type = PyArray_SHORT }; };

    template<> struct c_to_py<int> 
    { enum { py_type = PyArray_INT }; };

    template<> struct c_to_py<long> 
    { enum { py_type = PyArray_LONG }; };

    template<> struct c_to_py<float> 
    { enum { py_type = PyArray_FLOAT }; };

    template<> struct c_to_py<double> 
    { enum { py_type = PyArray_DOUBLE }; };

    template<> struct c_to_py<std::complex<float> >
    { enum { py_type = PyArray_CFLOAT }; };

    template<> struct c_to_py<std::complex<double> >
    { enum { py_type = PyArray_CDOUBLE }; };

    template<> struct c_to_py<PyObject*>
    { enum { py_type = PyArray_OBJECT }; };


    // class to get C type corresponding to PyArray_XXX enum type constant

    template<int py_type> struct py_to_c {};

    template<> struct py_to_c<PyArray_CHAR> 
    { typedef char c_type; };

    template<> struct py_to_c<PyArray_UBYTE> 
    { typedef unsigned char c_type; };

    template<> struct py_to_c<PyArray_SBYTE> 
    { typedef signed char c_type; };

    template<> struct py_to_c<PyArray_SHORT> 
    { typedef short c_type; };

    template<> struct py_to_c<PyArray_INT> 
    { typedef int c_type; };

    template<> struct py_to_c<PyArray_LONG> 
    { typedef long c_type; };

    template<> struct py_to_c<PyArray_FLOAT> 
    { typedef float c_type; };

    template<> struct py_to_c<PyArray_DOUBLE> 
    { typedef double c_type; };

    template<> struct py_to_c<PyArray_CFLOAT> 
    { typedef std::complex<float> c_type; };

    template<> struct py_to_c<PyArray_CDOUBLE> 
    { typedef std::complex<double> c_type; };

    template<> struct py_to_c<PyArray_OBJECT> 
    { typedef PyObject* c_type; };

    // get string with the PyArray_XXX mnemonic corresponding
    // to C type - for use in error reporting code

    template<class c_type> const char* c_to_py_str() {
      return type_info::get_s_type(c_to_py<c_type>::py_type);
    }

  } // namespace pyarray

class array : public object
{

public:

  array();

  explicit array(ref p);

  // construct the new array based on:
  // int ndim - number of dimensions
  // int dims[ndim] - number of elements in each dimension
  // T_elem type_el - sets the type of elements in the array, actual value does not matter

  template<class T_elem> explicit array(std::size_t ndim, const int *dims, T_elem type_el):
  object(ref(PyArray_FromDims(static_cast<int>(ndim), 
			      const_cast<int*>(dims), 
			      pyarray::c_to_py<T_elem>::py_type))) {}

  // Version of a constructor that uses external data rather than allocates its own storage.
  // This is of cause unsafe: the externally supplied data should not be deallocated
  // during lifetime of this object. NumPy manual claimes that it means it should be
  // never deallocated because "lifetime of Python objects is difficult to predict".

  template<class T_elem> explicit array(std::size_t ndim, const int *dims, T_elem *_data):
  object(ref(PyArray_FromDimsAndData(static_cast<int>(ndim), 
				     const_cast<int*>(dims), 
				     pyarray::c_to_py<T_elem>::py_type,
				     reinterpret_cast<char*>(_data)))) {}


  static PyTypeObject* type_obj();

  static bool accepts(ref p);

  // make a copy of this object with it's own data

  array clone() const;

  // wrapper around PyArray_Take() function -
  // see NumPy docs for the full description, it is rather complex. 
  // If, for example, this array has rank 2 (a matrix),
  // and 'indices' is 1-d array of 'int's,
  // the function will return matrix composed of rows
  // whose indices are listed in 'indices' if axes==0 or
  // composed from columns selected correspondingly if axes==1
  // Note: the NumPy docs dont' say it, but the look at the Numeric source
  // shows that PyArray_Take() always returns new non-shared copy of the 
  // input data. Therefore, this function cannot be used to select view
  // at the 'region of interest' in original array.
  // This function accepts indices as 'array' since PyArray_Take()
  // converts it's PyObject *indices argument to PyArray_Object anyway.

  array take(array indices, int axis);


  // return PyArray_XXX element type stored in this array

  int species() const;

  // rank of this array

  int rank() const;

  // number of elements in i-th dimension, index 'i' is zero-based

  int dimension(int i) const;

  // pointer to array with dimensions

  const int* dimensions() const;

  // writes dimensions to the output iterator

  template<typename _outiter> _outiter get_dimensions(_outiter res) const {
    const int *adim = dimensions();
    for(int i = 0; i < rank(); ++i, ++res)
      *res = adim[i];
    return res;
  }

  // 'strides' of this array as stored inside NumPy array instance,
  // i.e. in units of size 'char'

  const int* strides_char() const;

  // write strides of this array to output iterator;
  // unlike in 'strides_char()', strides are expressed in
  // units of the data elements stored in this array, not 'char'.
  // This method is expected to be used in conjuction with 'data_cast()'
  // method to get both the pointer to array of particular C data type
  // (e.g. 'float*') and strides for that datatype.

  template<class _OutIter> 
  _OutIter get_strides(_OutIter el_strides) const {
    int type_num = data_type_num();
    std::size_t el_size = pyarray::type_info::get_size(type_num);
    std::size_t size_div = el_size/sizeof(char);
    // Since we take data type from NumPy array itself,
    // we don't have to worry that size_div is a round devider
    // of char_strides[i].
    const int *char_strides = strides_char();
    for(int i = 0; i < rank(); ++i, ++el_strides)
      *el_strides = char_strides[i]/size_div;
    return el_strides;
  }

  // another way to get the PyArray_XXX type constant (see species())
  // (TODO: eliminate one of these functions)

  int data_type_num() const;

  // get string with mnemonic name of type constant for data in this array 
  // e.g "PyArray_FLOAT"

  const char* data_type_num_str() const {
    return pyarray::type_info::get_s_type(data_type_num());
  }

  // number of elements

  std::size_t size() const;

  bool is_contiguous() const;

  // pointer to char* data array stored in this NumPy object

  char* data_char() const;

  // cast this array's data pointer to a desired type,
  // throw an exception if requested data type does not
  // correspond to type_num value of this array (PyArray_XXX)

  // 'T_el dummy' argument is used only to specify the type to convert to,
  // the value of 'dummy' does not matter

  template<typename T_el> T_el* data_cast(T_el dummy) const {

    const int requested_py_num = pyarray::c_to_py<T_el>::py_type;
    const int data_py_num = data_type_num();

    if( ! is_data_type(dummy) ) {

      // If 'int' has same size as 'long', we should not distinguish
      // between NumPy arrays with PyArray_LONG and PyArray_INT data.
    if( ! ((
	    (requested_py_num == PyArray_LONG && data_py_num == PyArray_INT) ||
	    (requested_py_num == PyArray_INT && data_py_num == PyArray_LONG)
	    ) && sizeof(long) == sizeof(int) ) 
	) {
      std::string errmsg("Expected NumPy array with elements of type ");
      errmsg += pyarray::c_to_py_str<T_el>();
      errmsg += ". Instead, received array with elements of type ";
      errmsg += data_type_num_str();
      errmsg += ".";
      PyErr_SetString(PyExc_TypeError, errmsg.c_str());
      throw error_already_set();
    }
    }
    
    return reinterpret_cast<T_el*>(data_char());

  }

  // return contiguously allocated version of this array;
  // if this array is already contiguous, return the reference
  // to this array, otherwise, create and return contiguous copy

  array as_contiguous() const;

  // check that C data type stored in this array is indeed T_el

  template<typename T_el> bool is_data_type(T_el) const {
    return pyarray::c_to_py<T_el>::py_type == data_type_num();
  }

};

  // the recommended species if converting to array

  int species(const object& ob);

  int species(PyObject *pyob);

  // A set of functions to convert Python objects to arrays.
  // Functions differ in that some always return contiguously allocated version
  // and others always return the array with it's own copy of data.
  // This distinctions are important when the argument is itself an array.
  // These functions are simple wrappers around several NumPy conversion functions,
  // with an effort to provide a systematic naming convention reflecting the properties
  // of their output. The implementations try to minimize the number of array allocations
  // nessesary to achieve the desired results.
  //
  // The "int min_dim = 1, int max_dim = 0" default arguments are used to specify the minimum
  // and maximum allowed number of dimensions in the resulting array. Value of '0' means no
  // restrictions.

  // ------------------------------------------------------------------------------------------
  // This subset of functions specifies the desired type of elements in returned array
  // as PyArray_XXX type constant (with PyArray_NOTYPE meaning the smallest acceptable type)

  // if 'p' is already an array of desired type, return new reference to it
  // else return new array made from 'p' with it's own copy of data

  array to_array(PyObject* p,int min_dim = 1, int max_dim = 0,int i_type = PyArray_NOTYPE);

  // same as above, but will 'downcast' the data type if nessesary disregarding the possible 
  // loss of accuracy (e.g. cast 'float' to 'int' )

  array to_array_downcast(PyObject* p,int min_dim = 1, int max_dim = 0,int i_type = PyArray_NOTYPE);

  // always return new contiguous array made from 'p' and with it's own copy of data

  array to_array_copy(PyObject* p,int min_dim = 1, int max_dim = 0,int i_type = PyArray_NOTYPE);

  // same as above, but will 'downcast' the data type if nessesary disregarding the possible 
  // loss of accuracy (e.g. cast 'float' to 'int' )

  array to_array_copy_downcast(PyObject* p,int min_dim = 1, int max_dim = 0,int i_type = PyArray_NOTYPE);

  // if 'p' is an array and is contiguous, return new reference to it,
  // else return new array with it's own data

  array to_array_contiguous(PyObject* p,int min_dim = 1, int max_dim = 0,int i_type = PyArray_NOTYPE);

  // same as above, but will 'downcast' the data type if nessesary disregarding the possible 
  // loss of accuracy (e.g. cast 'float' to 'int' )

  array to_array_contiguous_downcast(PyObject* p,int min_dim = 1, int max_dim = 0,int i_type = PyArray_NOTYPE);

  //---------------------------------------------------------------------------------------------------

  // This subset of functions will take template parameter to determine data type of elements in returned array,
  // with the same semantics as corresponding non-template functions above

  template<typename T_el> 
  array to_array(PyObject *p,int min_dim = 1,int max_dim = 0) {  
    return to_array(p,min_dim,max_dim,pyarray::c_to_py<T_el>::py_type);
  }

  template<typename T_el> 
  array to_array_downcast(PyObject *p,int min_dim = 1,int max_dim = 0) {  
    return to_array_downcast(p,min_dim,max_dim,pyarray::c_to_py<T_el>::py_type);
  }

  template<typename T_el> 
  array to_array_copy(PyObject *p,int min_dim = 1,int max_dim = 0) {  
    return to_array_copy(p,min_dim,max_dim,pyarray::c_to_py<T_el>::py_type);
  }

  template<typename T_el> 
  array to_array_copy_downcast(PyObject *p,int min_dim = 1,int max_dim = 0) {  
    return to_array_copy_downcast(p,min_dim,max_dim,pyarray::c_to_py<T_el>::py_type);
  }

  template<typename T_el> 
  array to_array_contiguous(PyObject *p,int min_dim = 1,int max_dim = 0) {  
    return to_array_contiguous(p,min_dim,max_dim,pyarray::c_to_py<T_el>::py_type);
  }

  template<typename T_el> 
  array to_array_contiguous_downcast(PyObject *p,int min_dim = 1,int max_dim = 0) {  
    return to_array_contiguous_downcast(p,min_dim,max_dim,pyarray::c_to_py<T_el>::py_type);
  }

  //---------------------------------------------------------------------------------------


  // convenience function - check that rank of
  // array is N and throw descriptive Python exception
  // if not.
  
  bool assert_rank(const array& pyarr, int N, const char* add_msg=0);

} } // namespace boost::python


BOOST_PYTHON_BEGIN_CONVERSION_NAMESPACE

PyObject* to_python(const boost::python::array&);
boost::python::array from_python(PyObject* p, boost::python::type<boost::python::array>);

inline boost::python::array from_python(PyObject* p, boost::python::type<const boost::python::array&>)
{
    return from_python(p, boost::python::type<boost::python::array>());
}

BOOST_PYTHON_END_CONVERSION_NAMESPACE

#endif // ARRAY_BOOST_PY_H_

