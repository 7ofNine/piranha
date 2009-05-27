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

#ifndef PIRANHA_Q_EXPO_ARRAY_H
#define PIRANHA_Q_EXPO_ARRAY_H

#include <boost/algorithm/string/split.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/tuple/tuple.hpp>
#include <cmath>
#include <memory> // For default allocator.
#include <string>
#include <vector>

#include "../base_classes/q_array.h"
#include "../base_classes/toolbox.h"
#include "../common_functors.h"
#include "../exceptions.h"
#include "../mp.h"
#include "../q_power_cache.h"
#include "../utils.h"

#define __PIRANHA_Q_EXPO_ARRAY_TP_DECL int Pos, class Allocator
#define __PIRANHA_Q_EXPO_ARRAY_TP Pos,Allocator

namespace piranha
{
	template < __PIRANHA_Q_EXPO_ARRAY_TP_DECL = std::allocator<char> >
	struct q_expo_array {
		typedef toolbox<q_expo_array<__PIRANHA_Q_EXPO_ARRAY_TP> > type;
	};

	/// Rational exponents array.
	/**
	 * It wraps a piranha::q_array and adds the
	 * capabilities needed for exponent manipulation.
	 */
	template < __PIRANHA_Q_EXPO_ARRAY_TP_DECL >
	class toolbox<q_expo_array<__PIRANHA_Q_EXPO_ARRAY_TP> >: public q_array<__PIRANHA_Q_EXPO_ARRAY_TP, toolbox<q_expo_array<__PIRANHA_Q_EXPO_ARRAY_TP> > >
	{
			typedef q_array<__PIRANHA_Q_EXPO_ARRAY_TP, toolbox<q_expo_array<__PIRANHA_Q_EXPO_ARRAY_TP> > > ancestor;
			friend class q_array<__PIRANHA_Q_EXPO_ARRAY_TP, toolbox<q_expo_array<__PIRANHA_Q_EXPO_ARRAY_TP> > >;
			template <class SubSeries, class ArgsTuple>
			class sub_cache: public q_power_cache<SubSeries, typename base_series_arithmetics<SubSeries,ArgsTuple>::type>
			{
					typedef q_power_cache<SubSeries,
						typename base_series_arithmetics<SubSeries,ArgsTuple>::type> ancestor;
				public:
					sub_cache():ancestor() {}
					void setup(const SubSeries &s, const ArgsTuple *args_tuple) {
						this->m_arith_functor.m_args_tuple = args_tuple;
						this->m_container[mp_rational(0)] = SubSeries().base_add(1,*args_tuple);
						this->m_container[mp_rational(1)] = s;
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
			toolbox(): ancestor() {}
			/// Ctor from string.
			template <class ArgsTuple>
			explicit toolbox(const std::string &s, const ArgsTuple &): ancestor() {
				std::vector<std::string> sd;
				boost::split(sd, s, boost::is_any_of(std::string(1, this->separator)));
				const size_type w = boost::numeric_cast<size_type>(sd.size());
				this->resize(w);
				for (size_t i = 0; i < w; ++i) {
					(*this)[i] = utils::lexical_converter<value_type>(sd[i]);
				}
			}
			/// Ctor from psym.
			template <class ArgsTuple>
			explicit toolbox(const psym &p, const int &n, const ArgsTuple &a): ancestor(p, n, a) {}
			/// Multiplication.
			template <class QExpoArray, class ResultType>
			void multiply(const QExpoArray &qe2, ResultType &ret) const
			{
				const size_type max_w = this->size(), min_w = qe2.size();
				// Resize, if needed.
				ret.resize(max_w);
				// Assert widths, *this should always come from a polynomial, and its width should hence be
				// already adjusted my merge_args in multiplication routines.
				piranha_assert(max_w >= min_w);
				piranha_assert(ret.size() == max_w);
				size_type i;
				for (i = 0; i < min_w; ++i) {
					ret[i] = (*this)[i];
					ret[i] += qe2[i];
				}
				for (; i < max_w; ++i) {
					ret[i] = (*this)[i];
				}
			}
			// I/O.
			template <class ArgsTuple>
			void print_plain(std::ostream &out_stream, const ArgsTuple &args_tuple) const
			{
				piranha_assert(args_tuple.template get<ancestor::position>().size() == this->size());
				(void)args_tuple;
				this->print_elements(out_stream);
			}
			template <class ArgsTuple>
			void print_pretty(std::ostream &out_stream, const ArgsTuple &args_tuple) const
			{
				piranha_assert(args_tuple.template get<ancestor::position>().size() == this->size());
				bool printed_something = false;
				for (size_t i = 0; i < this->size(); ++i) {
					const value_type &q = (*this)[i];
					// Don't print anything if q is zero.
					if (q != 0) {
						// Prepend the multiplication operator only if we already printed something.
						if (printed_something) {
							out_stream << '*';
						}
						// Take care of printing the name of the exponent.
						out_stream << args_tuple.template get<ancestor::position>()[i].get_name();
						// Print the pow operator only if exponent is not unitary.
						if (q != 1) {
							out_stream << "**";
							const bool need_bracket = (q.get_den() != 1);
							if (need_bracket) {
								out_stream << "(";
							}
							out_stream << q;
							if (need_bracket) {
								out_stream << ")";
							}
						}
						printed_something = true;
					}
				}
			}
			template <class ArgsTuple>
			double eval(const double &t, const ArgsTuple &args_tuple) const
			{
				const size_t w = this->size();
				piranha_assert(w <= args_tuple.template get<ancestor::position>().size());
				double retval = 1.;
				for (size_t i = 0; i < w; ++i) {
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
			bool is_ignorable(const ArgsTuple &) const
			{
				return false;
			}
			bool is_unity() const
			{
				return (this->elements_are_zero());
			}
			/// Equality test.
			bool operator==(const toolbox &qe2) const
			{
				return this->elements_equal_to(qe2);
			}
			/// Norm.
			/**
			 * The norm of an exponent array is defined as the absolute value of the evaluation at t=0.
			 */
			template <class ArgsTuple>
			double norm(const ArgsTuple &args_tuple) const
			{
				return std::abs(eval(0, args_tuple));
			}
			/// Calculate hash value.
			size_t hash_value() const
			{
				return this->elements_hasher();
			}
			/// Return the total degree of the exponents array.
			int degree() const
			{
				value_type tmp(0);
				for (size_type i = 0; i < this->size(); ++i) {
					tmp += (*this)[i];
				}
				return sup_integer(tmp);
			}
			/// Total degree of the variables at specified positions pos.
			/**
			 * pos_tuple must be a tuple of vectors of (bool,size_t) pairs.
			 */
			template <class PosTuple>
			int partial_degree(const PosTuple &pos_tuple) const
			{
				const std::vector<std::pair<bool,size_t> > &pos = pos_tuple.template get<ancestor::position>();
				const size_type w = this->size(), pos_size = boost::numeric_cast<size_type>(pos.size());
				value_type tmp(0);
				for (size_type i = 0; i < pos_size; ++i) {
					// Add up exponents only if they are present and don't try to read outside boundaries
					// (this last could happen after merging arguments with a second series with smaller
					// number of arguments).
					if (pos[i].first && pos[i].second < w) {
						tmp += (*this)[pos[i].second];
					}
				}
				return sup_integer(tmp);
			}
			/// Minimum total degree.
			/**
			 * Provided for use within the power series toolbox, and defined to be equivalent to degree().
			 */
			int min_degree() const
			{
				return degree();
			}
			/// Minimum total degree of the variables at specified positions pos.
			/**
			 * Provided for use within the power series toolbox, and defined to be equivalent to partial_degree().
			 */
			template <class PosTuple>
			int partial_min_degree(const PosTuple &pos_tuple) const
			{
				return partial_degree(pos_tuple);
			}
			/// Return the position of the linear argument in the monomial.
			/**
			 * It will throw if the monomial is not linear or zero degree.
			 */
			int linear_arg_position() const
			{
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
			Series partial(const PosTuple &pos_tuple, const ArgsTuple &args_tuple) const
			{
				piranha_assert(pos_tuple.template get<ancestor::position>().size() == 1);
				const size_t pos = pos_tuple.template get<ancestor::position>()[0].second;
				piranha_assert(pos < this->size());
				// Do something only if the argument of the partial derivation is present in the exponent array
				// and the interesting exponent is not zero.
				Series retval;
				if (pos_tuple.template get<ancestor::position>()[0].first && (*this)[pos] != 0) {
					toolbox copy(*this);
					// TODO: use decrement operator here.
					copy[pos] = copy[pos] - 1;
					retval = Series::base_series_from_key(copy,args_tuple);
					retval.base_mult_by((*this)[pos],args_tuple);
				}
				return retval;
			}
			/// Exponentiation to double.
			template <class ArgsTuple>
			toolbox pow(const double &x, const ArgsTuple &args_tuple) const
			{
				return pow_impl(x,args_tuple);
			}
			/// Exponentiation to rational.
			template <class ArgsTuple>
			toolbox pow(const mp_rational &q, const ArgsTuple &args_tuple) const
			{
				return pow_impl(q,args_tuple);
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
					const size_t pos = pos_tuple.template get<ancestor::position>()[0].second;
					piranha_assert(pos < this->size());
					toolbox tmp_ea(*this);
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
			static int sup_integer(const mp_rational &q)
			{
				mp_integer num(q.get_num()), den(q.get_den());
				if (den == 1) {
					// If degree is already an integer number, let's
					// convert num to int and return it.
					return num.to_int();
				} else {
					// Divide truncating the non-integer part.
					num /= den;
					const int retval = num.to_int();
					// Take the ceil.
					return (q > 0) ? retval + 1 : retval;
				}
			}
			template <class T, class ArgsTuple>
			toolbox pow_impl(const T &x, const ArgsTuple &) const
			{
				toolbox retval(*this);
				if (x == -1) {
					retval.invert_sign();
				} else {
					const size_type w = this->size();
					for (size_type i = 0; i < w; ++i) {
						retval[i] *= x;
					}
				}
				return retval;
			}
	};
}

#undef __PIRANHA_Q_EXPO_ARRAY_TP_DECL
#undef __PIRANHA_Q_EXPO_ARRAY_TP

#endif
