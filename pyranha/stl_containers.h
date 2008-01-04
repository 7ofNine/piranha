/***************************************************************************
 *   Copyright (C) 2007, 2008 by Francesco Biscani   *
 *   bluescarni@gmail.com   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#ifndef PYRANHA_STL_CONTAINERS_H
#define PYRANHA_STL_CONTAINERS_H

#include <boost/python/class.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/return_value_policy.hpp>
#include <exception>

/// Functions to wrap STL vectors as Python lists.
template<class T>
  struct vector_to_list_helpers
{
  typedef typename T::value_type V;
  static V get(const T &x, int i)
  {
      if(i < 0) i+=x.size();
      if(i >= 0 and i < (int)x.size()) return x[i];
      IndexError();
      return 0;
  }
  static void set(T &x, int i, const V &v)
  {
      if(i < 0) i+=x.size();
      if(i >= 0 and i < (int)x.size()) x[i]=v;
      else IndexError();
  }
  static void del(T &x, int i)
  {
      if(i < 0) i+=x.size();
      if(i >= 0 and i < (int)x.size()) x.erase(x.begin()+i);
      else IndexError();
  }
  static void add(T &x, const V &v)
  {
      x.push_back(v);
  }
// This exception will be automatically translated by boost.python.
  static void IndexError() {throw(std::out_of_range("List index out of range."));}
};

/// Wrap STL vector as read write Python list.
template <class T>
  inline void vector_to_list(const std::string &name, const std::string &descr)
{
  boost::python::class_<T> inst(name.c_str(),descr.c_str());
  inst.def("__len__", &T::size);
  inst.def("clear", &T::clear);
  inst.def("append", &vector_to_list_helpers<T>::add,
    boost::python::with_custodian_and_ward<1,2>()); // to let container keep value
  inst.def("__getitem__", &vector_to_list_helpers<T>::get);
  inst.def("__setitem__", &vector_to_list_helpers<T>::set,
    boost::python::with_custodian_and_ward<1,3>()); // to let container keep value
  inst.def("__delitem__", &vector_to_list_helpers<T>::del);
}

/// Wrap STL vector as read-only Python list.
template <class T>
  inline void vector_to_rolist(const std::string &name, const std::string &descr)
{
  boost::python::class_<T> inst(name.c_str(),descr.c_str());
  inst.def("__len__", &T::size);
  inst.def("__getitem__", &vector_to_list_helpers<T>::get);
}

#endif
