/***************************************************************************
 *   Copyright (C) 2007 by Francesco Biscani   *
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

#include <boost/algorithm/string/trim.hpp>
#include <boost/lexical_cast.hpp>
#include <iostream>
#include <fstream>

#include "common_typedefs.h"
#include "settings_manager.h"
#include "stream_manager.h"

namespace piranha
{
  class utils
  {
    public:
/// Lexical converter.
/**
 * Convert a string to type T using boost::lexical_cast. If the operation is unsuccessful
 * a defaul constructed value is returned.
 * @param[in] s std::string to be converted.
 */
      template <class T>
        static T lexical_converter(const std::string &s)
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
      static int get_valid_string(std::ifstream &inf, std::string &str)
      {
        do
        {
          if (getline(inf,str,'\n').eof())
          {
            return 1;
          }
          boost::trim(str);
        }
        while (is_invalid_string(str));
        return 0;
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
        boost::split(split_v,str,boost::is_any_of(stream_manager::data_separator()));
        vector_double retval(split_v.size());
        for (unsigned int j=0;j<split_v.size();++j)
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
        for (size_t i=0;i<v.size();++i)
        {
          retval.append(boost::lexical_cast<std::string>(v[i]));
          if (i!=v.size()-1)
          {
            retval.append(stream_manager::data_separator());
          }
        }
        return retval;
      }
    private:
/// Check whether a string is empty or commented.
      static bool is_invalid_string(std::string &str)
      {
        if (str=="")
        {
          return true;
        }
        if (str[0]=='#')
        {
          return true;
        }
        return false;
      }
  };

}
#endif
