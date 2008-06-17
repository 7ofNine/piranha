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

#ifndef PIRANHA_EXPO_ARRAY_COMMONS_H
#define PIRANHA_EXPO_ARRAY_COMMONS_H

#include <boost/algorithm/string/split.hpp>
#include <boost/static_assert.hpp>
#include <cmath> // For std::pow, most likely temporary.
#include <iostream>
#include <string>
#include <vector>

#include "../exceptions.h"
#include "../integer_typedefs.h"
#include "../psym.h"
#include "../settings.h"

#define derived_const_cast (static_cast<Derived const *>(this))
#define derived_cast (static_cast<Derived *>(this))

namespace piranha
{
	/// Common class for dense exponent array.
	/**
	 * Intended to extend piranha::int_array for the manipulation of exponent
	 * parts in polynomials.
	 */
	template <class Derived>
	class expo_array_commons
	{
		public:
			// I/O.
			template <class ArgsTuple>
			void print_plain(std::ostream &out_stream, const ArgsTuple &args_tuple) const {
				// We assert like this because we want to make sure we don't go out of boundaries,
				// and because in case of fixed-width we may have smaller size of v wrt to "real" size.
				p_assert(args_tuple.template get<Derived::position>().size() <= derived_const_cast->size());
				(void)args_tuple;
				derived_const_cast->print_elements(out_stream);
			}
			void print_latex(std::ostream &out_stream, const vector_psym_p &v) const {
				// TODO: implement.
			}
			template <class ArgsTuple>
			double eval(const double &t, const ArgsTuple &args_tuple) const {
				const size_t w = args_tuple.template get<Derived::position>().size();
				p_assert(w <= derived_const_cast->size());
				double retval = 1.;
				for (size_t i = 0;i < w;++i) {
					retval *= std::pow(args_tuple.template get<Derived::position>()[i]->eval(t), (*derived_const_cast)[i]);
				}
				return retval;
			}
			/// Am I ignorable?
			/**
			 * By construction an array of exponents is never ignorable.
			 */
			template <class ArgsTuple>
			bool is_ignorable(const ArgsTuple &) const {
				return false;
			}
			bool is_unity() const {
				return (derived_const_cast->elements_are_zero());
			}
			/// Equality test.
			bool operator==(const Derived &e2) const {
				return derived_const_cast->elements_equal_to(e2);
			}
			/// Less than.
			bool operator<(const Derived &e2) const {
				return derived_const_cast->elements_lex_comparison(e2);
			}
			/// Norm.
			/**
			 * The norm of an exponent array is always one.
			 */
			template <class ArgsTuple>
			double norm(const ArgsTuple &args_tuple) const {
				p_assert(args_tuple.template get<Derived::position>().size() <= derived_const_cast->size());
				(void)args_tuple;
				return 1;
			}
			/// Calculate hash value.
			size_t hash_value() const {
				return derived_const_cast->elements_hasher();
			}
			/// Return the total degree of the exponents array.
			max_fast_int degree() const {
				max_fast_int retval = 0;
				for (typename Derived::size_type i = 0; i < derived_const_cast->m_size; ++i) {
					retval += (*derived_const_cast)[i];
				}
				return retval;
			}
			// This is the min total degree of a collection.
			// In this case the collection has a single element, hence the minimum degree is the degree itself.
			max_fast_int min_degree() const {
				return degree();
			}
			/// Return the position of the linear argument in the monomial.
			/**
			 * It will throw if the monomial is not linear or zero degree.
			 */
			max_fast_int linear_arg_position() const {
				bool found_linear = false;
				bool is_unity = true;
				max_fast_int candidate = -1;
				for (typename Derived::size_type i = 0; i < derived_const_cast->m_size; ++i) {
					if ((*derived_const_cast)[i] == 1) {
						is_unity = false;
						if (found_linear) {
							found_linear = false;
							break;
						} else {
							candidate = i;
							found_linear = true;
						}
					} else if ((*derived_const_cast)[i] != 0) {
						is_unity = false;
						found_linear = false;
						break;
					}
				}
				if (!found_linear && !is_unity) {
					throw unsuitable("Monomial is not linear.");
				}
				return candidate;
			}
			/// Calculate partial derivative.
			/**
			 * Result is a pair consisting of an integer and an exponent array.
			 */
			template <class PosTuple, class ArgsTuple>
			std::pair<int, Derived> partial(const PosTuple &pos_tuple, const ArgsTuple &) const {
				std::pair<int, Derived> retval(0, Derived());
				const size_t pos = pos_tuple.template get<Derived::position>().second;
				p_assert(pos < derived_const_cast->size());
				// Do something only if the argument of the partial derivation is present in the exponent array
				// and the interesting exponent is not zero.
				// Otherwise the above retval will return, and it will deliver a zero integer multiplier to be multiplied
				// by the coefficient in the partial derivation of the whole term.
				if (pos_tuple.template get<Derived::position>().first && derived_const_cast->m_ptr[pos] != 0) {
					retval.second = *derived_const_cast;
					retval.first = derived_const_cast->m_ptr[pos];
					--retval.second[pos];
				}
				return retval;
			}
			/// Integer exponentiation.
			/**
			 * If the exponent array cannot be raised to the desired power, an exception will be thrown.
			 */
			template <class ArgsTuple>
			Derived pow(const max_fast_int &n, const ArgsTuple &) const {
				typedef typename Derived::size_type size_type;
				Derived retval(*derived_const_cast);
				const size_type w = derived_const_cast->size();
				// Integer power. Retval has already been set to this, modify integers in-place.
				for (size_type i = 0; i < w; ++i) {
					retval[i] *= (typename Derived::value_type)n;
				}
				return retval;
			}
			/// Real exponentiation.
			/**
			 * If the exponent array cannot be raised to the desired power, an exception will be thrown.
			 */
			template <class ArgsTuple>
			Derived pow(const double &, const ArgsTuple &) const {
				// Real power is ok only if expo_array is unity.
				if (!is_unity()) {
					throw unsuitable("Cannot raise non-unity exponent array to real power.");
				}
				return Derived(*derived_const_cast);
			}
			template <class ArgsTuple>
			Derived root(const max_fast_int &n, const ArgsTuple &args_tuple) const {
				if (n == 0) {
					throw division_by_zero();
				} else if (n == 1) {
					return Derived(*derived_const_cast);
				}
				return pow(1. / (double)(n), args_tuple);
			}
			void upload_min_exponents(std::vector<max_fast_int> &v) const {
				derived_const_cast->upload_ints_to(v);
			}
			void test_min_exponents(std::vector <max_fast_int> &v) const {
				derived_const_cast->test_min_ints(v);
			}
			// Return true if the exponents are smaller than those specified in the limits vector.
			template <class ArgsTuple>
			bool test_expo_limits(const std::vector<std::pair<size_t, max_fast_int> > &v, const ArgsTuple &) const {
				const size_t size = v.size();
				for (size_t i = 0; i < size; ++i) {
					p_assert(v[i].first < derived_const_cast->m_size);
					switch ((*derived_const_cast)[v[i].first] > v[i].second) {
					case true:
						return false;
					case false:
						;
					}
				}
				return true;
			}
		protected:
			expo_array_commons() {}
			explicit expo_array_commons(const std::string &s) {
				typedef typename Derived::value_type value_type;
				std::vector<std::string> sd;
				boost::split(sd, s, boost::is_any_of(std::string(1, derived_const_cast->separator)));
				// TODO: check here that we are not loading too many multipliers, outside expo_size_t range.
				// TODO: do it everywhere!
				const size_t w = sd.size();
				derived_cast->resize(w);
				for (size_t i = 0;i < w;++i) {
					(*derived_cast)[i] = utils::lexical_converter<value_type>(sd[i]);
				}
			}
	};
}

#undef derived_const_cast
#undef derived_cast

#endif
