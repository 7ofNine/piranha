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

#ifndef PIRANHA_UTILS_H
#define PIRANHA_UTILS_H

#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/lexical_cast.hpp>
#include <iostream>
#include <fstream>
#include <valarray>

#include "common_typedefs.h"
#include "settings_manager.h"
#include "stream_manager.h"

namespace piranha
{
  /// Transform from Class2 to Class1.
  /**
   * Non-specialized version, it will create a copy of the converted class.
   */
  template <class Class2, class Class1> class class_converter
  {
    public:
      /// Constructor.
      /**
       * @param[in] c class to be converted.
       */
      explicit class_converter(const Class2 &c) :
      result(c)
      {
      }
      /// Copy of the converted class.
      const Class1 result;
  };

  /// Specialized class converter.
  /**
   * It will be invoked when the type to convert from is the same as the converted type. A reference
   * to the convertee is stored inside the class.
   */
  template <class Class2> class class_converter<Class2, Class2>
  {
    public:
      /// Constructor.
      /**
       * @param[in] c class to be converted.
       */
      explicit class_converter(const Class2 &c) :
      result(c)
      {
      }
      /// Reference to the converted class.
      const Class2 &result;
  };

  template <bool AssignZero, class VectorType> struct layout_assign_helper
  {
    static void run(VectorType &v1, const VectorType &v2, const size_t &i)
    {
      p_assert(i < v2.size());
      v1[i]=v2[i];
    }
  };

  template <class VectorType> struct layout_assign_helper<true, VectorType>
  {
    static void run(VectorType &v1, const VectorType &, const size_t &i)
    {
      v1[i]=0;
    }
  };

  class utils
  {
    public:
      /// Lexical converter.
      /**
       * Convert a string to type T using boost::lexical_cast. If the operation is unsuccessful
       * a defaul constructed value is returned.
       * @param[in] s std::string to be converted.
       */
      template <class T> static T lexical_converter(const std::string &s)
      {
        T retval;
        try
        {
          retval=boost::lexical_cast<T>(s);
        }
        catch(boost::bad_lexical_cast &)
        {
          std::cout << "Error in lexical_converter, returning default ctor." << std::endl;
          retval=T();
        }
        return retval;
      }
      /// Read a valid string from a stream.
      /**
       * If the string is empty or a comment (i.e., it starts with '#'), try to fetch another string
       * from stream. The returned string has leading and trailing spaces removed.
       * @param[in] inf input stream.
       * @param[out] str std::string which will contain output string.
       */
      static bool get_valid_string(std::ifstream &inf, std::string &str)
      {
        do
        {
          if (inf.eof())
          {
            return false;
          }
          getline(inf, str, '\n');
          boost::trim(str);
        } while (!is_valid(str));
        return true;
      }
      /// Open a file searching also in the 'theories of motion' directory.
      /**
       * This function returns an empty string if the file was not found, a string with the complete
       * file path if the file was found.
       * @param[in] fn std::string filename to be opened.
       * @param[out] std::ifstream file will be opened on.
       */
      static std::string open_file(const std::string &fn, std::ifstream &inf)
      {
        std::string filename=fn;
        inf.open(filename.c_str());
        if (inf.fail())
        {
          inf.close();
          // Clear ifstream's state (it is done automatically on close() on Linux, but not on
          // Windows).
          inf.clear();
          filename=(settings_manager::theories_path()+std::string("/")+filename);
          inf.open(filename.c_str());
          if (inf.fail())
          {
            std::cout << "Error opening file \"" << fn << "\"." << std::endl;
            std::cout << "Tried also \"" << filename << "\"." << std::endl;
            inf.close();
            return std::string();
          }
          std::cout << "Found \"" << fn << "\" in theories of motion." << std::endl;
        }
        return filename;
      }
      /// Convert a string into a vector of numerical values.
      /**
       * @param[in] str std::string to be converted.
       */
      static vector_double str_to_vector_double(const std::string &str)
      {
        deque_string split_v;
        boost::split(split_v, str, boost::is_any_of(stream_manager::data_separator()));
        vector_double retval(split_v.size());
        for (unsigned int j=0; j<split_v.size(); ++j)
        {
          retval[j]=lexical_converter<double>(split_v[j]);
        }
        return retval;
      }
      /// Convert a vector of numerical values into a string
      /**
       * @param[in] v vector_double to be converted.
       */
      static std::string vector_double_to_str(const vector_double &v)
      {
        std::string retval;
        for (size_t i=0; i<v.size(); ++i)
        {
          retval.append(boost::lexical_cast<std::string>(v[i]));
          if (i!=v.size()-1)
          {
            retval.append(stream_manager::data_separator());
          }
        }
        return retval;
      }
      /// Cache all pointers to elements of a container into an array.
      template <class Container, class Contained>
        static void array_pointer(const Container &c, std::valarray<Contained const *> &v)
      {
        typedef typename Container::const_iterator const_iterator;
        const size_t l=c.length();
        v.resize(l);
        size_t i=0;
        const const_iterator it_f=c.end();
        for (const_iterator it=c.begin();it!=it_f;++it)
        {
          v[i]=&(*it);
          ++i;
        }
      }
      /// Get relative layout of two vectors.
      template <class VectorType>
        static layout_type get_layout(const VectorType &v1, const VectorType &v2)
      {
        const size_t size1 = v1.size(), size2 = v2.size();
        // First we must construct v2's layout wrt to v1.
        layout_type l(size2);
        for (size_t i=0;i < size2;++i)
        {
          // If we won't find v2's element, we'll mark it as not found.
          l[i].first=false;
          // For each of ps2's elements, look for that same element in v1.
          for (size_t j=0;j < size1;++j)
          {
            if (v1[j] == v2[i])
            {
              // We found it, mark as found and proceed to next v2 element.
              l[i].first=true;
              l[i].second=j;
              break;
            }
          }
        }
        // Now we must take care of those elements of v1 that are not represented in the layout (i.e., they are not in v2)
        for (size_t i=0;i < size1;++i)
        {
          // Look for element index i in the layout.
          bool found = false;
          const size_t l_size = l.size();
          for (size_t j=0;j < l_size;++j)
          {
            if (l[j].first and l[j].second == i)
            {
              found = true;
              break;
            }
          }
          // If we did not find it, append it to the layout.
          if (!found)
          {
            l.push_back(layout_element(true,i));
          }
        }
        return l;
      }
      /// Apply layout to vector of arguments.
      /**
       * Applies layout l to vector v1, where the missing elements in the layout are taken from v2.
       */
      template <class VectorType>
        static void apply_layout(const layout_type &l, VectorType &v1, const VectorType &v2)
      {
        generic_apply_layout<false,VectorType>(l,v1,v2);
      }
      /// Apply layout to vector of arguments.
      /**
       * Applies layout l to vector v1, where the missing elements are zeroed.
       */
      template <class VectorType>
        static void apply_layout(const layout_type &l, VectorType &v1)
      {
        generic_apply_layout<true,VectorType>(l,v1,VectorType());
      }
    private:
      /// Check whether a string is valid.
      /**
       * Invalid strings are empty or commented.
       */
      static bool is_valid(const std::string &str)
      {
        // TODO: replace with switch statement.
        if (str.empty() or str[0] == '#')
        {
          return false;
        }
        return true;
      }
      template <bool AssignZero, class VectorType>
        static void generic_apply_layout(const layout_type &l, VectorType &v1, const VectorType &v2)
      {
        const size_t l_size = l.size();
        // The layout must have at least all arguments in v1.
        p_assert(l_size >= v1.size());
        // Memorize the old vector.
        const VectorType old(v1);
        // Make space.
        v1.resize(l_size);
        for (size_t i=0;i < l_size;++i)
        {
          switch (l[i].first)
          {
            case true:
              p_assert(l[i].second < old.size());
              v1[i]=old[l[i].second];
              break;
            case false:
              layout_assign_helper<AssignZero,VectorType>::run(v1,v2,i);
          }
        }
      }
  };
}
#endif
