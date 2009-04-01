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
#include <boost/numeric/conversion/cast.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "settings.h"

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
			template <class T> static T lexical_converter(const std::string &s) {
				T retval;
				try {
					retval = boost::lexical_cast<T>(s);
				} catch (boost::bad_lexical_cast &) {
					std::cout << "Error in lexical_converter, returning default-constructed object." << '\n';
					retval = T();
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
			static bool get_valid_string(std::ifstream &inf, std::string &str) {
				do {
					if (inf.eof()) {
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
			static std::string open_file(const std::string &fn, std::ifstream &inf) {
				std::string filename = fn;
				inf.open(filename.c_str(), std::ios::in | std::ios::binary);
				if (inf.fail()) {
					inf.close();
					// Clear ifstream's state (it is done automatically on close() on Linux, but not on
					// Windows).
					inf.clear();
					filename = (settings::path() + std::string("/") + filename);
					inf.open(filename.c_str(), std::ios::in | std::ios::binary);
					if (inf.fail()) {
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
			template <class T>
			static std::vector<T> str_to_vector(const std::string &str, const std::string &separator) {
				if (str.empty()) {
					return std::vector<T>();
				}
				std::vector<std::string> v;
				boost::split(v, str, boost::is_any_of(separator));
				const size_t size = v.size();
				std::vector<T> retval(size);
				for (size_t j = 0; j < size; ++j) {
					retval[j] = lexical_converter<T>(v[j]);
				}
				return retval;
			}
			/// Convert a vector of numerical values into a string
			/**
			 * @param[in] v vector_double to be converted.
			 */
			template <class T>
			static std::string vector_to_str(const std::vector<T> &v, const std::string &separator) {
				std::string retval;
				const size_t size = v.size();
				for (size_t i = 0; i < size; ++i) {
					retval.append(boost::lexical_cast<std::string>(v[i]));
					if (i != size - 1) {
						retval.append(separator);
					}
				}
				return retval;
			}
			/// Check if a double-precision number is an integer.
			static bool is_integer(const double &x) {
				try {
					const int tmp = boost::numeric_cast<int>(x);
					if (x == tmp) {
						return true;
					}
				} catch (const boost::bad_numeric_cast &) {
					// This means that x is overflowing int's range.
				}
				return false;
			}
		private:
			/// Check whether a string is valid.
			/**
			 * Invalid strings are empty or commented.
			 */
			static bool is_valid(const std::string &str) {
				if (str.empty() || str[0] == '#') {
					return false;
				}
				return true;
			}
	};
}
#endif
