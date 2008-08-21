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

#ifndef PIRANHA_PSYMBOL_H
#define PIRANHA_PSYMBOL_H

#include <boost/array.hpp>
#include <boost/functional/hash.hpp>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/mem_fun.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <cmath>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "config.h"
#include "exceptions.h"
#include "p_assert.h"
#include "settings.h"
#include "utils.h"

namespace piranha
{
	class __PIRANHA_VISIBLE psyms
	{
		public:
			/// Literal symbol class.
			/**
			  * This class is used represent symbolic arguments. It features a string representing the symbol's name and a
			  * numerical vector
			  * which is used to evaluate the symbol in time in a polynomial fashion. For instance, if the numerical vector
			  * has a size of three and its elements are named \f$ \alpha \f$, \f$ \beta \f$ and \f$ \gamma \f$,
			  * it means that the symbol is evaluated as
			  * \f[
			  * \alpha + \beta t + \gamma t^2,
			  * \f]
			  * where \f$ t \f$ is time.
			  */
			class __PIRANHA_VISIBLE psym
			{
					friend class psyms;
				public:
					// Ctors
					/// Constructor from two std::string objects
					/**
					  * Assigns both psym::m_name and psym::m_time_eval by parsing the two strings.
					  * @param[in] name symbol's name.
					  * @param[in] te symbol's time evaluation expressed in string format.
					  */
					psym(const std::string &name, const std::string &te): m_name(name),
							m_time_eval(utils::str_to_vector<double>(te, separator)) {
						reg(*this);
					}
					/// Constructor from std::string.
					/**
					* Constructs a psym with empty time evaluation
					*/
					psym(const std::string &str): m_name(str), m_time_eval() {
						reg(*this);
					}
					// Constructors from multiple values for time evlauation.
					psym(const std::string &s, const double &x1): m_name(s), m_time_eval((size_t)1) {
						boost::array<double, 1> tmp = {
							{
								x1
							}
						};
						build_from_array(tmp);
					}
					psym(const std::string &s, const double &x1, const double &x2): m_name(s), m_time_eval((size_t)2) {
						boost::array<double, 2> tmp = {
							{
								x1, x2
							}
						};
						build_from_array(tmp);
					}
					psym(const std::string &s, const double &x1, const double &x2, const double &x3):
							m_name(s), m_time_eval((size_t)3) {
						boost::array<double, 3> tmp = {
							{
								x1, x2, x3
							}
						};
						build_from_array(tmp);
					}
					psym(const std::string &s, const double &x1, const double &x2, const double &x3, const double &x4):
							m_name(s), m_time_eval((size_t)4) {
						boost::array<double, 4> tmp = {
							{
								x1, x2, x3, x4
							}
						};
						build_from_array(tmp);
					}
					psym(const std::string &s, const double &x1, const double &x2, const double &x3,
						 const double &x4, const double &x5): m_name(s), m_time_eval((size_t)5) {
						boost::array<double, 5> tmp = {
							{
								x1, x2, x3, x4, x5
							}
						};
						build_from_array(tmp);
					}
					/// Print to stream.
					void print(std::ostream &out_stream = std::cout) const {
						settings::setup_stream(out_stream);
						out_stream << "name=" << m_name << '\n';
						out_stream << "time_eval=";
						for (size_t j = 0; j < m_time_eval.size(); ++j) {
							out_stream << m_time_eval[j];
							if (j != m_time_eval.size() - 1) {
								out_stream << separator;
							}
						}
						out_stream << '\n';
					}
					/// Time evaluation.
					double eval(const double &t) const {
						double retval = 0.;
						const size_t w = m_time_eval.size();
						for (size_t i = 0; i < w; ++i) {
							// FIXME: use natural_pow or the like here, to speed up?
							retval += std::pow(t, (int)i) * m_time_eval[i];
						}
						return retval;
					}
					std::string name() const {
						return m_name;
					}
					const std::vector<double> &time_eval() const {
						return m_time_eval;
					}
				private:
					// Helper for ctor from boost::array.
					template <class T>
					void build_from_array(const T &a) {
						const size_t size = a.size();
						p_assert(a.size() == m_time_eval.size());
						for (size_t i = 0; i < size; ++i) {
							m_time_eval[i] = a[i];
						}
						reg(*this);
					}
				private:
					// Data members.
					const std::string				m_name;
					// Mutable because we want to be able to freely change it in the psym manager.
					mutable std::vector<double>		m_time_eval;
					static const std::string		separator;
			};
		private:
			typedef boost::multi_index_container
			<
			psym,
			boost::multi_index::indexed_by
			<
			boost::multi_index::ordered_unique<boost::multi_index::const_mem_fun<psym, std::string, &psym::name> >
			>
			> set_type;
		public:
			typedef set_type::iterator psym_p;
			typedef psym_p iterator;
			static psym_p get_pointer(const psym &psym) {
				psym_p retval(get_pointer(psym.m_name));
				return retval;
			}
			static psym_p get_pointer(const std::string &name) {
				psym_p retval(set.find(name));
				if (retval == set.end()) {
					throw not_existing(std::string("Symbol \"") + name + "\" does not exist.");
				}
				return retval;
			}
			/// Print list of symbols.
			static void print(std::ostream &stream = std::cout) {
				const psym_p it_f = set.end();
				for (psym_p it = set.begin(); it != it_f; ++it) {
					it->print(stream);
				}
			}
			static psym get(const std::string &name) {
				return *get_pointer(name);
			}
			// Needed to mimic standard container.
			static size_t length() {
				return set.size();
			}
			static iterator begin() {
				return set.begin();
			}
			static iterator end() {
				return set.end();
			}
		private:
			static void reg(const psym &psym) {
				const psym_p it = set.find(psym.m_name);
				if (it == set.end()) {
					// Symbol is not already present, add it.
					std::pair<psym_p, bool> result = set.insert(psym);
					p_assert(result.second);
				} else {
					// Symbol name is already registered, overwrite time evaluation.
					it->m_time_eval = psym.m_time_eval;
				}
			}
		private:
			static set_type set;
	};

	/// Typedefs used in series, terms, coefficients and trigonometric parts.
	typedef psyms::psym psym;
	typedef psyms::psym_p psym_p;
	typedef std::vector<psym_p> vector_psym_p;
}
#endif
