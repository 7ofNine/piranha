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

#ifndef PIRANHA_EXPO_ARRAY_H
#define PIRANHA_EXPO_ARRAY_H

#include <boost/algorithm/string/split.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/numeric/conversion/converter.hpp>
#include <boost/type_traits/integral_constant.hpp> // For type traits.
#include <cmath> // For std::abs and std::pow (this last is most likely temporary).
#include <cstddef>
#include <iostream>
#include <memory> // For standard allocator.
#include <string>
#include <utility> // For std::pair.
#include <vector>

#include "../base_classes/int_array.h"
#include "../base_classes/vector_key.h"
#include "../common_functors.h"
#include "../exceptions.h"
#include "../int_power_cache.h"
#include "../math.h"
#include "../mp.h"
#include "../psym.h"
#include "../settings.h"
#include "../type_traits.h"
#include "../utils.h" // For lexical_cast.

#define __PIRANHA_EXPO_ARRAY_TP_DECL int Bits, int Pos, class Allocator
#define __PIRANHA_EXPO_ARRAY_TP Bits,Pos,Allocator

namespace piranha
{
	/// Exponents array.
	/**
	 * It wraps a piranha::int_array with integer sized Bits, and adds the
	 * capabilities needed for exponent manipulation.
	 */
	template < __PIRANHA_EXPO_ARRAY_TP_DECL = std::allocator<char> >
	class expo_array: public int_array<__PIRANHA_EXPO_ARRAY_TP, expo_array<__PIRANHA_EXPO_ARRAY_TP> >
	{
			typedef int_array<__PIRANHA_EXPO_ARRAY_TP, expo_array<__PIRANHA_EXPO_ARRAY_TP> > ancestor;
			friend class int_array<__PIRANHA_EXPO_ARRAY_TP, expo_array<__PIRANHA_EXPO_ARRAY_TP> >;
			template <class SubSeries, class ArgsTuple>
			class sub_cache: public int_power_cache<SubSeries, base_series_arithmetics<SubSeries,ArgsTuple> >
			{
					typedef int_power_cache<SubSeries,
						base_series_arithmetics<SubSeries,ArgsTuple> > ancestor;
				public:
					sub_cache():ancestor() {}
					void setup(const SubSeries &s, const ArgsTuple *args_tuple) {
						this->m_arith_functor.m_args_tuple = args_tuple;
						this->m_container[0] = SubSeries().base_add(1,*args_tuple);
						this->m_container[1] = s;
					}
			};
			// This is just a stub.
			template <class SubSeries, class ArgsTuple>
			class ei_sub_cache
			{
				public:
					void setup(const SubSeries &, const ArgsTuple *) {}
			};
		public:
			typedef int degree_type;
			typedef typename ancestor::value_type value_type;
			typedef typename ancestor::size_type size_type;
			typedef double eval_type;
			template <class SubSeries, class SubCachesCons, class ArgsTuple>
			struct sub_cache_selector {
				typedef boost::tuples::cons<sub_cache<SubSeries,ArgsTuple>,
					SubCachesCons> type;
			};
			template <class SubSeries, class SubCachesCons, class ArgsTuple>
			struct ei_sub_cache_selector {
				typedef boost::tuples::cons<ei_sub_cache<SubSeries,ArgsTuple>,
					SubCachesCons> type;
			};
			// Ctors.
			/// Default ctor.
			expo_array(): ancestor() {}
			/// Ctor from string.
			template <class ArgsTuple>
			explicit expo_array(const std::string &s, const ArgsTuple &): ancestor() {
				std::vector<std::string> sd;
				boost::split(sd, s, boost::is_any_of(std::string(1, this->separator)));
				// TODO: check here that we are not loading too many multipliers, outside expo_size_t range.
				// TODO: do it everywhere!
				const std::size_t w = sd.size();
				this->resize(w);
				for (std::size_t i = 0; i < w; ++i) {
					(*this)[i] = utils::lexical_converter<value_type>(sd[i]);
				}
			}
			/// Ctor from psym.
			template <class ArgsTuple>
			explicit expo_array(const psym &p, const int &n, const ArgsTuple &a): ancestor(p, n, a) {}
			// Math.
			/// Multiplication.
			template <class ExpoArray, class ResultType>
			void multiply(const ExpoArray &e2, ResultType &ret) const {
				const size_type max_w = this->size(), min_w = e2.size();
				// Resize, if needed.
				ret.resize(max_w);
				// Assert widths, *this should always come from a polynomial, and its width should hence be
				// already adjusted my merge_args in multiplication routines.
				piranha_assert(max_w >= min_w);
				piranha_assert(ret.size() == max_w);
				size_type i;
				for (i = 0;i < min_w;++i) {
					ret[i] = (*this)[i] + e2[i];
				}
				for (;i < max_w;++i) {
					ret[i] = (*this)[i];
				}
			}
			// I/O.
			template <class ArgsTuple>
			void print_plain(std::ostream &out_stream, const ArgsTuple &args_tuple) const {
				piranha_assert(args_tuple.template get<ancestor::position>().size() == this->size());
				(void)args_tuple;
				this->print_elements(out_stream);
			}
			template <class ArgsTuple>
			void print_pretty(std::ostream &out_stream, const ArgsTuple &args_tuple) const {
				piranha_assert(args_tuple.template get<ancestor::position>().size() == this->size());
				bool printed_something = false;
				for (std::size_t i = 0; i < this->m_size; ++i) {
					const int n = (*this)[i];
					// Don't print anything if n is zero.
					if (n != 0) {
						// Prepend the multiplication operator only if we already printed something.
						if (printed_something) {
							out_stream << '*';
						}
						// Take care of printing the name of the exponent.
						out_stream << args_tuple.template get<ancestor::position>()[i].get_name();
						// Print the pow operator only if exponent is not unitary.
						if (n != 1) {
							out_stream << "**" << n;
						}
						printed_something = true;
					}
				}
			}
			template <class ArgsTuple>
			void print_tex(std::ostream &out_stream, const ArgsTuple &args_tuple) const {
				piranha_assert(args_tuple.template get<ancestor::position>().size() == this->size());
				for (std::size_t i = 0; i < this->m_size; ++i) {
					const int n = (*this)[i];
					// Don't print anything if n is zero.
					if (n != 0) {
						// Take care of printing the name of the exponent.
						out_stream << ' ' << args_tuple.template get<ancestor::position>()[i].get_name() << ' ';
						// Print the pow operator only if exponent is not unitary.
						if (n != 1) {
							out_stream << "^{" << n << '}';
						}
					}
				}
			}
			template <class ArgsTuple>
			double eval(const double &t, const ArgsTuple &args_tuple) const {
				const std::size_t w = this->size();
				piranha_assert(w <= args_tuple.template get<ancestor::position>().size());
				double retval = 1.;
				for (std::size_t i = 0;i < w;++i) {
					retval *= std::pow(args_tuple.template get<ancestor::position>()[i].eval(t),
						(*this)[i]);
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
				return (this->elements_are_zero());
			}
			/// Equality test.
			bool operator==(const expo_array &e2) const {
				return this->elements_equal_to(e2);
			}
			/// Inequality test.
			bool operator<(const expo_array &e2) const {
				return this->lex_comparison(e2);
			}
			/// Norm.
			/**
			 * The norm of an exponent array is defined as the absolute value of the evaluation at t=0.
			 */
			template <class ArgsTuple>
			double norm(const ArgsTuple &args_tuple) const {
				return std::abs(eval(0, args_tuple));
			}
			/// Calculate hash value.
			std::size_t hash_value() const {
				return this->elements_hasher();
			}
			/// Return the total degree of the exponents array.
			int degree() const {
				int retval = 0;
				for (size_type i = 0; i < this->m_size; ++i) {
					retval += (*this)[i];
				}
				return retval;
			}
			/// Total degree of the variables at specified positions pos.
			/**
			 * pos_tuple must be a tuple of vectors of (bool,std::size_t) pairs.
			 */
			template <class PosTuple>
			int partial_degree(const PosTuple &pos_tuple) const {
				const std::vector<std::pair<bool,std::size_t> > &pos = pos_tuple.template get<ancestor::position>();
				const size_type w = this->size(), pos_size = boost::numeric_cast<size_type>(pos.size());
				int retval = 0;
				for (size_type i = 0; i < pos_size; ++i) {
					// Add up exponents only if they are present and don't try to read outside boundaries
					// (this last could happen after merging arguments with a second series with smaller
					// number of arguments).
					if (pos[i].first && pos[i].second < w) {
						retval += (*this)[pos[i].second];
					}
				}
				return retval;
			}
			/// Minimum total degree of the exponents array.
			/**
			 * Provided for use within the power series toolbox, and defined to be equivalent to degree().
			 */
			int order() const {
				return degree();
			}
			/// Minimum total degree of the variables at specified positions pos.
			/**
			 * Provided for use within the power series toolbox, and defined to be equivalent to partial_degree().
			 */
			template <class PosTuple>
			int partial_order(const PosTuple &pos_tuple) const {
				return partial_degree(pos_tuple);
			}
			/// Return the position of the linear argument in the monomial.
			/**
			 * It will throw if the monomial is not linear or zero degree.
			 */
			int linear_arg_position() const {
				bool found_linear = false;
				bool is_unity = true;
				int candidate = -1;
				for (size_type i = 0; i < this->m_size; ++i) {
					if ((*this)[i] == 1) {
						is_unity = false;
						if (found_linear) {
							found_linear = false;
							break;
						} else {
							candidate = i;
							found_linear = true;
						}
					} else if ((*this)[i] != 0) {
						is_unity = false;
						found_linear = false;
						break;
					}
				}
				if (!found_linear && !is_unity) {
					piranha_throw(value_error,"monomial is not linear");
				}
				return candidate;
			}
			/// Calculate partial derivative.
			template <class Series, class PosTuple, class ArgsTuple>
			Series partial(const PosTuple &pos_tuple, const ArgsTuple &args_tuple) const {
				piranha_assert(pos_tuple.template get<ancestor::position>().size() == 1);
				const std::size_t pos = pos_tuple.template get<ancestor::position>()[0].second;
				piranha_assert(!pos_tuple.template get<ancestor::position>()[0].first || pos < this->size());
				// Do something only if the argument of the partial derivation is present in the exponent array
				// and the interesting exponent is not zero.
				Series retval;
				if (pos_tuple.template get<ancestor::position>()[0].first && (*this)[pos] != 0) {
					expo_array copy(*this);
					--copy[pos];
					retval = Series::base_series_from_key(copy,args_tuple);
					retval.base_mult_by((*this)[pos],args_tuple);
				}
				return retval;
			}
			template <class ArgsTuple>
			expo_array pow(const double &y, const ArgsTuple &) const {
				if (is_integer(y)) {
					return pow_int((int)y);
				} else {
					return pow_double(y);
				}
			}
			template <class ArgsTuple>
			expo_array pow(const mp_rational &q, const ArgsTuple &) const {
				expo_array retval(*this);
				const size_type size = this->size();
				for(size_type i = 0; i < size; ++i) {
					mp_rational tmp(q);
					tmp *= (*this)[i];
					if (tmp.get_den() != 1) {
						piranha_throw(value_error,"exponent is not suitable for the calculation of rational power");
					}
					retval[i] = boost::numeric_cast<value_type>(tmp.to_int());
				}
				return retval;
			}
			template <class RetSeries, class PosTuple, class SubCaches, class ArgsTuple>
			RetSeries sub(const PosTuple &pos_tuple, SubCaches &sub_caches,
				const ArgsTuple &args_tuple) const
			{
				RetSeries retval;
				// If the argument is not present here, the return series will have one term consisting
				// of a unitary coefficient and this very expo_array.
				// NOTE: for now we can substitute one symbol at a time.
				piranha_assert(pos_tuple.template get<ancestor::position>().size() == 1);
				if (!pos_tuple.template get<ancestor::position>()[0].first) {
					retval = RetSeries::base_series_from_key(*this, args_tuple);
				} else {
					const std::size_t pos = pos_tuple.template get<ancestor::position>()[0].second;
					piranha_assert(pos < this->size());
					expo_array tmp_ea(*this);
					// Let's turn off the exponent associated to the symbol we are substituting.
					tmp_ea[pos] = 0;
					RetSeries orig(RetSeries::base_series_from_key(tmp_ea, args_tuple));
					piranha_assert(retval.empty());
					// NOTICE: series multadd here?
					retval.base_add(orig, args_tuple);
					retval.base_mult_by(sub_caches.template get<ancestor::position>()
						[(*this)[pos]],
					args_tuple);
				}
				return retval;
			}
			template <class RetSeries, class PosTuple, class SubCaches, class ArgsTuple>
			RetSeries ei_sub(const PosTuple &, SubCaches &,
				const ArgsTuple &args_tuple) const {
				return RetSeries::base_series_from_key(*this, args_tuple);
			}
		private:
			/// Integer exponentiation.
			/**
			 * If the exponent array cannot be raised to the desired power, an exception will be thrown.
			 */
			expo_array pow_int(const int &n) const {
				expo_array retval(*this);
				const size_type w = this->size();
				// Integer power. Retval has already been set to this, modify integers in-place.
				// TODO: check for overflow?
				for (size_type i = 0; i < w; ++i) {
					retval[i] *= (value_type)n;
				}
				return retval;
			}
			/// Real exponentiation.
			/**
			 * If the exponent array cannot be raised to the desired power, an exception will be thrown.
			 */
			expo_array pow_double(const double &) const {
				// Real power is ok only if expo_array is unity.
				if (!is_unity()) {
					piranha_throw(value_error,"cannot raise non-unity exponent array to real power");
				}
				return expo_array(*this);
			}
	};

	/// is_ring_exact type trait specialisation for expo_array.
	template <__PIRANHA_EXPO_ARRAY_TP_DECL>
	struct is_ring_exact<expo_array<__PIRANHA_EXPO_ARRAY_TP> >: boost::true_type {};
}

#undef __PIRANHA_EXPO_ARRAY_TP_DECL
#undef __PIRANHA_EXPO_ARRAY_TP

#endif
