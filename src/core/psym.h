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
	// the global "dictionary" of Psym(s)
	struct PIRANHA_VISIBLE PsymManager {

		//implementation of a Psym stored in the PsymManager
		struct PIRANHA_VISIBLE PsymImpl {

			PsymImpl(const std::string &name, const std::vector<double> &time_eval = std::vector<double>()):
				m_name(name), m_time_eval(time_eval) {}

			// Print to stream.
			void print(std::ostream &outStream) const
			{
				outStream << "name=" << m_name << '\n';
				outStream << "time_eval=";
				for (std::size_t j = 0; j < m_time_eval.size(); ++j) 
				{
					outStream << boost::lexical_cast<std::string>(m_time_eval[j]);
					if (j != m_time_eval.size() - 1) 
					{
						outStream << separator;
					}
				}
				outStream << '\n';
			}


			bool operator<(const PsymImpl &other) const
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


			// Data members of PsymImpl
			const std::string		    m_name;
			// Mutable because we want to be able to freely change it in the Psym manager.
			mutable std::vector<double>	m_time_eval;
			static const std::string	separator;
		};

		// PsymManager global storage variable
		typedef std::set<PsymImpl> container_type;
		static container_type       container;
	};

	// Forward declaration for use in the typedef below.
	class Psym;

	typedef std::vector<Psym> VectorPsym;

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
	class PIRANHA_VISIBLE Psym {

			typedef PsymManager::container_type::const_iterator it_type;
			typedef PsymManager::PsymImpl PsymImpl;

			struct push_back_to 
			{
				push_back_to(VectorPsym &v):m_v(&v) {}


				template <class T>
				void operator()(const T &x) const
				{
					m_v->push_back(Psym(x.m_name, x.m_time_eval));
				}


				mutable VectorPsym *m_v;
			};

		public:

			/// Constructor from name and time evaluation in string form.
			Psym(const std::string &name, const std::string &time_eval)
			{
				const PsymImpl p(name, utils::str_to_vector<double>(time_eval, PsymImpl::separator));
				connstructFromImpl(p);
			}


			/// Constructor from name.
			/**
			 * If the symbol is already present in the Psym manager then use it, otherwise initialise
			 * a new Psym with the given name and an empty time evaluation vector.
			 */
			explicit Psym(const std::string &name)
			{
				PsymImpl const p(name);

				it_type const it = PsymManager::container.find(p);

				if (it == PsymManager::container.end()) 
				{
					std::pair<it_type, bool> const res = PsymManager::container.insert(p);
					PIRANHA_ASSERT(res.second);
					m_it = res.first;

				} else 
				{
					m_it = it;
				}
			}


			/// Constructor from name and vector.
			Psym(const std::string &name, const std::vector<double> &time_eval)
			{
				const PsymImpl p(name, time_eval);
				connstructFromImpl(p);
			}


			/// Constructor from name and constant value.
			Psym(const std::string &name, const double &value)
			{
				const PsymImpl p(name, std::vector<double>(std::size_t(1), value));
				connstructFromImpl(p);
			}


			bool operator<(const Psym &other) const
			{
				return ((*m_it) < *(other.m_it));
			}


			bool operator==(const Psym &other) const
			{
				return m_it == other.m_it;
			}


			bool operator!=(const Psym &other) const
			{
				return !(*this == other);
			}


			static VectorPsym list()
			{
				VectorPsym retval;
				retval.reserve(PsymManager::container.size());

				std::for_each(PsymManager::container.begin(), PsymManager::container.end(), push_back_to(retval));
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
// TODO:     this does not work with the Python interfaces. How to fix
//            void set_time_eval(const std::string &t) const
//            {
//                
//                m_it->m_time_eval = utils::str_to_vector<double>(t, PsymImpl::separator);
//            }

		private:

			void connstructFromImpl(const PsymImpl &p)
			{
				const it_type it = PsymManager::container.find(p);
				if (it == PsymManager::container.end()) 
				{
					const std::pair<it_type, bool> res = PsymManager::container.insert(p);
					PIRANHA_ASSERT(res.second);
					m_it = res.first;
				} else 
				{
					it->m_time_eval = p.m_time_eval;
					m_it = it;
				}
			}

		private:
			// a Psym is an iterator into the set of PsymImpl managed in a global set in PsymManager
			it_type m_it;
	};

	template <class ArgsTuple>
	struct psyms2pos_impl {
		template <class PosTuple>
		static void run(const VectorPsym &v, PosTuple &pos_tuple, const ArgsTuple &argsTuple)
		{
			const std::size_t a_size = argsTuple.get_head().size(), v_size = v.size();
			pos_tuple.get_head().reserve(v_size);
			// For each psymbol, test presence.
			for (std::size_t i = 0; i < v_size; ++i) 
			{
				// Initially set the symbol to not found.
				pos_tuple.get_head().push_back(std::make_pair(false,std::size_t(0)));
				for (std::size_t j = 0; j < a_size; ++j) 
				{
					if (argsTuple.get_head()[j] == v[i]) 
					{
						pos_tuple.get_head().back().first  = true;
						pos_tuple.get_head().back().second = j;
						// No need to continue, the symbol is (supposed to be) unique.
						break;
					}
				}
			}
			psyms2pos_impl<typename ArgsTuple::tail_type>::run(v,pos_tuple.get_tail(),argsTuple.get_tail());
		}
	};


	template <>
	struct psyms2pos_impl<boost::tuples::null_type> 
	{
		static void run(const VectorPsym &, const boost::tuples::null_type &, const boost::tuples::null_type &) {}
	};


	// Transform a vector of psyms into a tuple of vectors of (position flag, position) pairs, given a reference arguments tuple.
	// Return value will be a tuple of vectors, each of size v.size(), containing (presence, position) pairs for the corresponding symbols
	// in v.
	template <class ArgsTuple>
	inline typename NTuple<std::vector<std::pair<bool, std::size_t> >, boost::tuples::length<ArgsTuple>::value>::Type  psyms2pos(VectorPsym const &v, ArgsTuple const &argsTuple)
	{
		typedef typename NTuple<std::vector<std::pair<bool, std::size_t> >, boost::tuples::length<ArgsTuple>::value>::Type retval_type;

		// First we want to make sure that the vector of symbols does not contain duplicate elements.
		std::set<Psym> const uniquesSet(v.begin(), v.end());
		VectorPsym     const uniquesVector(uniquesSet.begin(), uniquesSet.end());
		retval_type retval;

		psyms2pos_impl<ArgsTuple>::run(uniquesVector, retval, argsTuple);

		return retval;
	}


	// Transform a vector of names into a vector of symbols.
	inline VectorPsym names2psyms(const std::vector<std::string> &vs)
	{
		const std::size_t size = vs.size();
		VectorPsym v;
		v.reserve(size);
		for (std::size_t i = 0; i < size; ++i) 
		{
			v.push_back(Psym(vs[i]));
		}
		return v;
	}
}
#endif
