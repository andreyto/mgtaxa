/*
   Copyright: Andrei Tovchigrechko, 2001 

   The implementation of wrapper class around NumPy Numeric array object for boost::python library.

*/


#include <boost/python/objects.hpp>
#include <boost/python/detail/none.hpp>

#include "DockTK/Pyarray/array.hpp"

#include <string>

#include "DockTK/Common/to_string.hpp"

namespace boost { namespace python {

  namespace pyarray {

    const type_info_rec type_info::atype_info[] = {
      {PyArray_CHAR,    sizeof(char),                 "PyArray_CHAR" },
      {PyArray_UBYTE,   sizeof(unsigned char),        "PyArray_UBYTE"},
      {PyArray_SBYTE,   sizeof(signed char),          "PyArray_SBYTE"},
      {PyArray_SHORT,   sizeof(short),                "PyArray_SHORT"},
      {PyArray_INT,     sizeof(int),                  "PyArray_INT"},
      {PyArray_LONG,    sizeof(long),                 "PyArray_LONG"},
      {PyArray_FLOAT,   sizeof(float),                "PyArray_FLOAT"},
      {PyArray_DOUBLE,  sizeof(double),               "PyArray_DOUBLE"},
      {PyArray_CFLOAT,  sizeof(std::complex<float>),  "PyArray_CFLOAT"},
      {PyArray_CDOUBLE, sizeof(std::complex<double>), "PyArray_CDOUBLE"},
      {PyArray_OBJECT,  sizeof(PyObject*),            "PyArray_OBJECT"}
    };

    std::size_t type_info::get_size(int itype) {

      for( int i = 0; i < PyArray_NTYPES; i++)
	if( atype_info[i].itype == itype ) return  atype_info[i].size;

      PyErr_SetString(PyExc_TypeError, "Unknown PyArray element type" );
      throw error_already_set();

      return 0;
    }

    const char* type_info::get_s_type(int itype) {

      for( int i = 0; i < PyArray_NTYPES; i++)
	if( atype_info[i].itype == itype ) return  atype_info[i].s_type;

      PyErr_SetString(PyExc_TypeError, "Unknown PyArray element type" );
      throw error_already_set();

      return 0;
    }

    struct TypeProperties {
      TypeProperties() {
	longIsInt = (type_info::get_size(PyArray_LONG) == type_info::get_size(PyArray_INT));
      }

      bool longIsInt;
    };

    TypeProperties typeProperties;

    ref fixType(ref p) {
      ref a(PyArray_FromObject(p.get(),PyArray_NOTYPE,0,0));
      if( PyArray_ObjectType( a.get(), 0) == PyArray_LONG && 
	  ( ! typeProperties.longIsInt ) ) {
	a = ref(PyArray_Cast((PyArrayObject*)a.get(),PyArray_INT));
      }
      return a;
    }

  } // namespace pyarray


  // hope that having object with null pointer inside does
  // not break things; existing objects in boost::python 
  // (e.g. dictionary) don't do it
  array::array() :
    object(ref((PyObject*)0,ref::null_ok))
  {}

array::array(ref p)
  : object(pyarray::fixType(p))
{

    assert(accepts(p));
    if (!accepts(p))
    {
        PyErr_SetString(PyExc_TypeError, p->ob_type->tp_name);
        throw error_already_set();
    }
}


PyTypeObject* array::type_obj()
{
    return &PyArray_Type;
}

bool array::accepts(ref p)
{
    return PyArray_Check(p.get());
}

int array::data_type_num() const {
  return ((PyArrayObject*) get())->descr->type_num;
}


array array::clone() const {
  PyObject *p = PyArray_CopyFromObject(get(), species(), rank(), rank());
  return array(ref(p));
}

array array::take(array indices, int axis) {
  return array(ref(PyArray_Take(get(),indices.get(),axis)));
}

int  array::species() const {
  return PyArray_ObjectType(get(), 0);
}

int  array::rank() const {
  return ((PyArrayObject*) get())->nd;
}

int  array::dimension(int i) const {
  assert(i>=0 && i<rank());
  return ((PyArrayObject*) get())->dimensions[i];
}

const int*  array::dimensions() const {
  return ((PyArrayObject*) get())->dimensions;
}

const int*  array::strides_char() const {
  return ((PyArrayObject*) get())->strides;
}

bool  array::is_contiguous() const {
  return PyArray_ISCONTIGUOUS ((PyArrayObject*) get()) == 1;
}

char*  array::data_char() const {
  return ((PyArrayObject*) get())->data;
}

array  array::as_contiguous() const {
  if (is_contiguous()) return array(reference());
  else return array(ref(
			PyArray_ContiguousFromObject(get(), species(), rank(), rank())
			));
}        

std::size_t array::size() const {
  return static_cast<std::size_t>(PyArray_Size(get()));
}


int species(PyObject *pyob) {
    return PyArray_ObjectType(pyob, 0);
}

array to_array(PyObject* p,int min_dim, int max_dim,int i_type) {
  return array(ref(PyArray_FromObject(p,i_type,min_dim,max_dim)));
}

array to_array_downcast(PyObject* p,int min_dim, int max_dim, int i_type) {
    PyObject *p_arr = PyArray_FromObject(p, i_type, min_dim, max_dim);
    if( p_arr == NULL ) {
      PyErr_Clear(); // clear the error condition created by PyArray_CopyFromObject()
      array arr = to_array(p,min_dim,max_dim,PyArray_NOTYPE);
      p_arr = PyArray_Cast((PyArrayObject*)arr.get(),i_type);
    }
    return array(ref(p_arr));
}

array to_array_copy(PyObject* p,int min_dim, int max_dim,int i_type) {
    return array(ref(PyArray_CopyFromObject(p,i_type, min_dim, max_dim)));
}

array to_array_copy_downcast(PyObject* p,int min_dim, int max_dim, int i_type) {
    PyObject *p_arr = PyArray_CopyFromObject(p, i_type, min_dim, max_dim);
    if( p_arr == NULL ) {
      PyErr_Clear(); // clear the error condition created by PyArray_CopyFromObject()
      // since first call failed, cast will certainly be performed if at all possible in PyArray_Cast,
      // hence it is safe to call just to_array() - PyArray_Cast will return new contiguous non-shared array
      array arr = to_array(p,min_dim,max_dim,PyArray_NOTYPE);
      p_arr = PyArray_Cast((PyArrayObject*)arr.get(),i_type);
    }
    return array(ref(p_arr));
}


array to_array_contiguous(PyObject* p,int min_dim, int max_dim, int i_type) {
    return array(ref(PyArray_ContiguousFromObject(p, i_type, min_dim, max_dim)));
}

array to_array_contiguous_downcast(PyObject* p,int min_dim, int max_dim, int i_type) {
    PyObject *p_arr = PyArray_ContiguousFromObject(p, i_type, min_dim, max_dim);
    if( p_arr == NULL ) {
      PyErr_Clear(); // clear the error condition created by PyArray_ContiguousFromObject()
      // since first call failed, cast will certainly be performed if at all possible in PyArray_Cast,
      // hence it is safe to call just to_array() - PyArray_Cast will return new contiguous non-shared array 
      array arr = to_array(p,min_dim,max_dim,PyArray_NOTYPE);
      p_arr = PyArray_Cast((PyArrayObject*)arr.get(),i_type);
    }
    return array(ref(p_arr));
}

bool assert_rank(const array& pyarr, int N, const char *add_msg) {
    // assert that rank NumPy array is N
    if( pyarr.rank() != N ) {

      std::string errmsg("Rank of supplied NumPy array is ");
      errmsg += DockTK::to_string(pyarr.rank());
      errmsg += ". Rank expected by our code is ";
      errmsg += DockTK::to_string(N);
      errmsg += ". ";
      if( add_msg ) errmsg += add_msg;
      PyErr_SetString(PyExc_TypeError,errmsg.c_str());
      throw error_already_set();
      return false;
    }
    return true;
}

}} // namespace boost::python


// this functions are defined in boost/libs/python/opjects.cpp
namespace boost { namespace python {

template <class T>
inline T object_from_python(PyObject* p, type<T>)
{
    ref x(p, ref::increment_count);
    if (!T::accepts(x))
    {
        PyErr_SetString(PyExc_TypeError, p->ob_type->tp_name);
        throw error_already_set();
    }
    return T(x);
}

inline PyObject* object_to_python(const object& x)
{
    return x.reference().release();
}

}} // namespace boost::python

BOOST_PYTHON_BEGIN_CONVERSION_NAMESPACE

PyObject* to_python(const boost::python::array& x)
{
    return boost::python::object_to_python(x);
}

boost::python::array from_python(PyObject* p, boost::python::type<boost::python::array> type)
{
    return boost::python::object_from_python(p, type);
}

BOOST_PYTHON_END_CONVERSION_NAMESPACE

