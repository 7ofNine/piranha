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

#include <algorithm>
#include <boost/lexical_cast.hpp>
#include <cstddef>
#include <iostream>
#include <set>
#include <string>
#include <utility> // For std::pair.
#include <vector>

#include <cmath>

#include "config.h"
#include "exceptions.h"
#include "ntuple.h"
#include "settings.h"
#include "utils.h"

namespace piranha
{
	struct __PIRANHA_VISIBLE psym_manager {
		struct __PIRANHA_VISIBLE psym_impl {

			psym_impl(const std::string &name, const std::vector<double> &time_eval = std::vector<double>()):
				m_name(name),m_time_eval(time_eval) {}

			// Print to stream.
			void print(std::ostream &out_stream) const
			{
				out_stream << "name=" << m_name << '\n';
				out_stream << "time_eval=";
				for (std::size_t j = 0; j < m_time_eval.size(); ++j) 
				{
					out_stream << boost::lexical_cast<std::string>(m_time_eval[j]);
					if (j != m_time_eval.size() - 1) 
					{
						out_stream << separator;
					}
				}
				out_stream << '\n';
			}


			bool operator<(const psym_impl &other) const
			{
				return (m_name < other.m_name);
			}


			// Time evaluation.
			double eval(const double &t) const
			{
				double retval = 0.;
				const std::size_t w = m_time_eval.size();
				for (std::size_t i = 0; i < w; ++i) 
				{
					// FIXME: use natural_pow or the like here, to speed up?
					retval += std::pow(t, (int)i) * m_time_eval[i];
				}
				return retval;
			}


			// Data members.
			const std::string		    m_name;
			// Mutable because we want to be able to freely change it in the psym manager.
			mutable std::vector<double>	m_time_eval;
			static const std::string	separator;
		};

		typedef std::set<psym_impl> container_type;
		static container_type       container;
	};

	// Forward declaration for use in the typedef below.
	class psym;

	typedef std::vector<psym> vector_psym;

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
	class __PIRANHA_VISIBLE psym {
			typedef psym_manager::container_type::const_iterator it_type;
			typedef psym_manager::psym_impl psym_impl;
			struct push_back_to {
				
				push_back_to(vector_psym &v):m_v(&v) {}

				template <class T>
				void operator()(const T &x) const
				{
					m_v->push_back(psym(x.m_name,x.m_time_eval));
				}

				mutable vector_psym *m_v;
			};

		public:

			/// Constructor from name and time evaluation in string form.
			explicit psym(const std::string &name, const std::string &time_eval)
			{
				const psym_impl p(name,utils::str_to_vector<double>(time_eval,psym_impl::separator));
				construct_from_impl(p);
			}


			/// Constructor from name.
			/**
			 * If the symbol is already present in the psym manager then use it, otherwise initialise
			 * a new psym with the given name and an empty time evaluation vector.
			 */
			explicit psym(const std::string &name)
			{
				const psym_impl p(name);
				const it_type it = psym_manager::container.find(p);
				if (it == psym_manager::container.end()) 
				{
					const std::pair<it_type,bool> res = psym_manager::container.insert(p);
					piranha_assert(res.second);
					m_it = res.first;
				} else 
				{
					m_it = it;
				}
			}


			/// Constructor from name and vector.
			explicit psym(const std::string &name, const std::vector<double> &time_eval)
			{
				const psym_impl p(name,time_eval);
				construct_from_impl(p);
			}


			/// Constructor from name and constant value.
			explicit psym(const std::string &name, const double &value)
			{
				const psym_impl p(name,std::vector<double>((std::size_t)1,value));
				construct_from_impl(p);
			}


			bool operator<(const psym &other) const
			{
				return ((*m_it) < *(other.m_it));
			}


			bool operator==(const psym &other) const
			{
				return m_it == other.m_it;
			}


			bool operator!=(const psym &other) const
			{
				return !(*this == other);
			}


			static vector_psym list()
			{
				vector_psym retval;
				retval.reserve(psym_manager::container.size());
				std::for_each(psym_manager::container.begin(),psym_manager::container.end(),push_back_to(retval));
				return retval;
			}


			// TODO: move to operator<< for streams? Along with other classes...
			/// Print to stream.
			void print(std::ostream &s = std::cout) const
			{
				m_it->print(s);
			}


			/// Evaluate symbol at time t.
			double eval(const double &t) const
			{
				return m_it->eval(t);
			}


			/// Name getter.
			const std::string &get_name() const
			{
				return m_it->m_name;
			}


			/// Time evaluation vector getter.
			const std::vector<double> &get_time_eval() const
			{
				return m_it->m_time_eval;
			}


			/// Time evaluation vector setter.
			void set_time_eval(const std::vector<double> &t) const
			{
				m_it->m_time_eval = t;
			}

		private:

			void construct_from_impl(const psym_impl &p)
			{
				const it_type it = psym_manager::container.find(p);
				if (it == psym_manager::container.end()) 
				{
					const std::pair<it_type,bool> res = psym_manager::container.insert(p);
					piranha_assert(res.second);
					m_it = res.first;
				} else 
				{
					it->m_time_eval = p.m_time_eval;
					m_it = it;
				}
			}

		private:
			
			it_type m_it;
	};

	template <class ArgsTuple>
	struct psyms2pos_impl {
		template <class PosTuple>
		static void run(const vector_psym &v, PosTuple &pos_tuple, const ArgsTuple &args_tuple)
		{
			const std::size_t a_size = args_tuple.get_head().size(), v_size = v.size();
			pos_tuple.get_head().reserve(v_size);
			// For each psymbol, test presence.
			for (std::size_t i = 0; i < v_size; ++i) 
			{
				// Initially set the symbol to not found.
				pos_tuple.get_head().push_back(std::make_pair(false,std::size_t(0)));
				for (std::size_t j = 0; j < a_size; ++j) 
				{
					if (args_tuple.get_head()[j] == v[i]) 
					{
						pos_tuple.get_head().back().first  = true;
						pos_tuple.get_head().back().second = j;
						// No need to continue, the symbol is (supposed to be) unique.
						break;
					}
				}
			}
			psyms2pos_impl<typename ArgsTuple::tail_type>::run(v,pos_tuple.get_tail(),args_tuple.get_tail());
		}
	};


	template <>
	struct psyms2pos_impl<boost::tuples::null_type> 
	{
		static void run(const vector_psym &, const boost::tuples::null_type &, const boost::tuples::null_type &) {}
	};


	// Transform a vector of psyms into a tuple of vectors of (position flag, position) pairs, given a reference arguments tuple.
	// Return value will be a tuple of vectors, each of size v.size(), containing (presence, position) pairs for the corresponding symbols
	// in v.
	template <class ArgsTuple>
	inline typename ntuple<std::vector<std::pair<bool,std::size_t> >,boost::tuples::length<ArgsTuple>::value>::type psyms2pos(const vector_psym &v, const ArgsTuple &args_tuple)
	{
		typedef typename ntuple<std::vector<std::pair<bool,std::size_t> >,boost::tuples::length<ArgsTuple>::value>::type retval_type;
		// First we want to make sure that the vector of symbols does not contain duplicate elements.
		const std::set<psym> uniques_set(v.begin(),v.end());
		const vector_psym uniques_vector(uniques_set.begin(),uniques_set.end());
		retval_type retval;
		psyms2pos_impl<ArgsTuple>::run(uniques_vector,retval,args_tuple);
		return retval;
	}


	// Transform a vector of names into a vector of symbols.
	inline vector_psym names2psyms(const std::vector<std::string> &vs)
	{
		const std::size_t size = vs.size();
		vector_psym v;
		v.reserve(size);
		for (std::size_t i = 0; i < size; ++i) 
		{
			v.push_back(psym(vs[i]));
		}
		return v;
	}
}
#endif
