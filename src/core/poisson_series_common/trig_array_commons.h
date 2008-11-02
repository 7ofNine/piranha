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

#ifndef PIRANHA_TRIG_ARRAY_COMMONS_H
#define PIRANHA_TRIG_ARRAY_COMMONS_H

#include <boost/algorithm/string/split.hpp>
#include <cmath> // For std::abs.
#include <complex>
#include <string>
#include <utility> // For std::pair.
#include <vector>

#include "../base_classes/series_builders.h"
#include "../config.h"
#include "../psym.h"
#include "../settings.h"
#include "trig_evaluator.h"

#define derived_const_cast (static_cast<Derived const *>(this))
#define derived_cast (static_cast<Derived *>(this))

namespace piranha
{
	/// Common class for dense trigonometric array.
	/**
	 * Intended to extend piranha::int_array for the manipulation of trigonometric
	 * parts in Poisson series.
	 */
	template <class Derived>
	class trig_array_commons
	{
		public:
			/// Return const reference to flavour.
			const bool &flavour() const {
				return derived_const_cast->m_flavour;
			}
			/// Return mutable reference to flavour.
			bool &flavour() {
				return derived_cast->m_flavour;
			}
			// I/O.
			template <class ArgsTuple>
			void print_plain(std::ostream &out_stream, const ArgsTuple &args_tuple) const {
				// We assert like this because we want to make sure we don't go out of boundaries,
				// and because in case of fixed-width we may have smaller size of v wrt to "real" size.
				p_assert(args_tuple.template get<Derived::position>().size() <= derived_const_cast->size());
				(void)args_tuple;
				derived_const_cast->print_elements(out_stream);
				// Print the separator before flavour only if we actually printed something above.
				if (derived_const_cast->size() != 0) {
					out_stream << derived_const_cast->separator;
				}
				if (derived_const_cast->m_flavour) {
					out_stream << "c";
				} else {
					out_stream << "s";
				}
			}
			template <class ArgsTuple>
			void print_pretty(std::ostream &out_stream, const ArgsTuple &args_tuple) const {
				p_assert(args_tuple.template get<Derived::position>().size() <= derived_const_cast->m_size);
				if (derived_const_cast->m_flavour) {
					out_stream << "cos(";
				} else {
					out_stream << "sin(";
				}
				bool printed_something = false;
				for (size_t i = 0; i < derived_const_cast->m_size; ++i) {
					const max_fast_int n = derived_const_cast->m_container.v[i];
					// Don't print anything if n is zero.
					if (n != 0) {
						// If we already printed something and n is positive we are going to print the sign too.
						if (printed_something && n > 0) {
							out_stream << '+';
						}
						// Take care of printing the multiplier.
						if (n == 1) {
							;
						} else if (n == -1) {
							out_stream << '-';
						} else {
							out_stream << n << '*';
						}
						out_stream << args_tuple.template get<Derived::position>()[i]->name();
						printed_something = true;
					}
				}
				out_stream << ')';
			}
			void print_latex(std::ostream &out_stream, const vector_psym_p &v) const {
				const size_t w = v.size();
				p_assert(w <= derived_const_cast->size())
				switch (derived_const_cast->m_flavour) {
				case true:
					out_stream << "c&";
					break;
				case false:
					out_stream << "s&";
				}
				bool first_one = true;
				std::string tmp("$");
				for (size_t i = 0;i < w;++i) {
					if ((*derived_const_cast)[i] != 0) {
						if ((*derived_const_cast)[i] > 0 && !first_one) {
							tmp.append("+");
						}
						if ((*derived_const_cast)[i] == -1) {
							tmp.append("-");
						} else if ((*derived_const_cast)[i] == 1) {} else {
							tmp.append(boost::lexical_cast<std::string>((int)(*derived_const_cast)[i]));
						}
						tmp.append(v[i]->name());
						first_one = false;
					}
				}
				tmp.append("$");
				// If we did not write anything erase math markers.
				if (tmp == "$$") {
					tmp.clear();
				}
				out_stream << tmp;
			}
			bool is_unity() const {
				return (derived_const_cast->elements_are_zero() && derived_const_cast->m_flavour);
			}
			/// Frequency.
			/**
			 * Get the frequency of the linear combination, given a vector of piranha::psym pointers describing
			 * the arguments.
			 * @param[in] v vector of piranha::psym pointers.
			 */
			template <class ArgsTuple>
			double freq(const ArgsTuple &args_tuple) const {
				return combined_time_eval<1>(args_tuple);
			}
			/// Phase.
			/**
			 * Get the phase of the linear combination, given a vector of piranha::psym pointers describing the
			 * arguments.
			 * @param[in] v vector of piranha::psym pointers.
			 */
			template <class ArgsTuple>
			double phase(const ArgsTuple &args_tuple) const {
				return combined_time_eval<0>(args_tuple);
			}
			/// Norm.
			/**
			 * The norm of a trigonometric part is always one.
			 */
			template <class ArgsTuple>
			double norm(const ArgsTuple &args_tuple) const {
				p_assert(args_tuple.template get<Derived::position>().size() >= derived_const_cast->size());
				(void)args_tuple;
				return 1;
			}
			/// Time evaluation of arguments.
			/**
			 * Returns the value assumed by the linear combination of arguments at time t.
			 * @param[in] t double time of the evaluation.
			 * @param[in] v vector of piranha::psym pointers.
			 */
			template <class ArgsTuple>
			double eval(const double &t, const ArgsTuple &args_tuple) const {
				const size_t w = derived_const_cast->size();
				p_assert(w <= args_tuple.template get<Derived::position>().size());
				double retval = 0.;
				for (size_t i = 0;i < w;++i) {
					if ((*derived_const_cast)[i] != 0) {
						retval += (*derived_const_cast)[i] * args_tuple.template get<Derived::position>()[i]->eval(t);
					}
				}
				if (derived_const_cast->m_flavour) {
					return std::cos(retval);
				} else {
					return std::sin(retval);
				}
			}
			/// Time evaluation of complex exponential of the arguments.
			/**
			 * Returns the real or imaginary part (depending on flavour) of the complex exponential of the
			 * linear combination of arguments at time t.
			 * Uses a piranha::trig_evaluator object which contains a cache of the complex exponentials of arguments.
			 * @param[in] te piranha::trig_evaluator containing a cache of complex exponentials of arguments.
			 */
			template <class TrigEvaluator>
			double t_eval(TrigEvaluator &te) const {
				const size_t w = te.width();
				p_assert(w <= derived_const_cast->size());
				std::complex<double> retval(1.);
				for (size_t i = 0;i < w;++i) {
					if ((*derived_const_cast)[i] != 0) {
						retval *= te.request_ei(i, (*derived_const_cast)[i]);
					}
				}
				if (derived_const_cast->m_flavour) {
					return retval.real();
				} else {
					return retval.imag();
				}
			}
			/// Sign.
			/**
			 * Return the sign of the first non-zero element of the combination. Zero is considered positive.
			 * This function is used to test for canonical form in piranha::poisson_series_term.
			 */
			short int sign() const {
				const size_t w = derived_const_cast->size();
				for (size_t i = 0;i < w;++i) {
					// TODO: use switch?
					if ((*derived_const_cast)[i] > 0) {
						return 1;
					}
					if ((*derived_const_cast)[i] < 0) {
						return -1;
					}
				}
				return 1;
			}
			// All multipliers are zero and flavour is sine.
			template <class ArgsTuple>
			bool is_ignorable(const ArgsTuple &) const {
				return (derived_const_cast->elements_are_zero() && !derived_const_cast->m_flavour);
			}
			/// Equality test.
			bool operator==(const Derived &t2) const {
				return (derived_const_cast->m_flavour == t2.m_flavour && derived_const_cast->elements_equal_to(t2));
			}
			/// Less than.
			bool operator<(const Derived &t2) const {
				if (derived_const_cast->m_flavour < t2.m_flavour) {
					return true;
				} else if (derived_const_cast->m_flavour > t2.m_flavour) {
					return false;
				}
				return derived_const_cast->lex_comparison(t2);
			}
			/// Calculate hash_value.
			/**
			 * Used by the hash_value overload for piranha::base_term.
			 */
			size_t hash_value() const {
				size_t retval = derived_const_cast->elements_hasher();
				boost::hash_combine(retval, derived_const_cast->m_flavour);
				return retval;
			}
			/// Multiply by an integer.
			void mult_by_int(const int &n) {
				const size_t w = derived_const_cast->size();
				for (size_t i = 0;i < w;++i) {
					(*derived_cast)[i] *= n;
				}
			}
			/// Calculate partial derivative.
			/**
			 * Result is a pair consisting of an integer and a trigonometric array.
			 */
			template <class PosTuple, class ArgsTuple>
			std::pair<int, Derived> partial(const PosTuple &pos_tuple, const ArgsTuple &) const {
				std::pair<int, Derived> retval(0, Derived());
				// Do something only if the argument of the partial derivation is present in the trigonometric array.
				// Otherwise the above retval will return, and it will deliver a zero integer multiplier to be
				// multiplied by the coefficient in the partial derivation of the whole term.
				if (pos_tuple.template get<Derived::position>().first) {
					retval.second = *derived_const_cast;
					const size_t pos = pos_tuple.template get<Derived::position>().second;
					// Change the flavour of the resulting key.
					retval.second.m_flavour = (!derived_const_cast->m_flavour);
					p_assert(pos < derived_const_cast->size());
					retval.first = derived_const_cast->m_container.v[pos];
					if (derived_const_cast->m_flavour) {
						retval.first = -retval.first;
					}
				}
				return retval;
			}
			template <class ArgsTuple>
			Derived pow(const max_fast_int &n, const ArgsTuple &) const {
				const bool int_zero = derived_const_cast->elements_are_zero();
				Derived retval;
				if (n < 0) {
					if (int_zero && !derived_const_cast->m_flavour) {
						// 0^-n.
						throw division_by_zero();
					} else if (int_zero && derived_const_cast->m_flavour) {
						// 1^-n == 1. Don't do nothing because retval is already initialized properly.
						;
					} else {
						// x^-n -> no go.
						throw unsuitable("Non-unity Trigonometric array is not suitable for negative integer "
							"exponentiation.");
					}
				} else if (n == 0) {
					// x^0 == 1. Don't do nothing because retval is already initialized properly.
					;
				} else {
					if (int_zero && !derived_const_cast->m_flavour) {
						// 0^n == 0.
						retval.m_flavour = false;
					} else if (int_zero && derived_const_cast->m_flavour) {
						// 1^y == 1. Don't do nothing because retval is already initialized properly.
						;
					} else {
						// x^n --> no go (it should be handled by natural power routine for series).
						throw unsuitable("Non-unity Trigonometric array is not suitable for positive integer"
							"exponentiation.");
					}
				}
				return retval;
			}
			/// Real exponentiation.
			template <class ArgsTuple>
			Derived pow(const double &y, const ArgsTuple &) const {
				const bool int_zero = derived_const_cast->elements_are_zero();
				Derived retval;
				if (y < 0) {
					if (int_zero && !derived_const_cast->m_flavour) {
						// 0^-y.
						throw division_by_zero();
					} else if (int_zero && derived_const_cast->m_flavour) {
						// 1^-y == 1. Don't do nothing because retval is already initialized properly.
						;
					} else {
						// x^-y -> no go.
						throw unsuitable("Non-unity Trigonometric array is not suitable for negative real "
							"exponentiation.");
					}
				} else if (y == 0) {
					// x^0 == 1. Don't do nothing because retval is already initialized properly.
					;
				} else {
					if (int_zero && !derived_const_cast->m_flavour) {
						// 0^y == 0.
						retval.m_flavour = false;
					} else if (int_zero && derived_const_cast->m_flavour) {
						// 1^y == 1. Don't do nothing because retval is already initialized properly.
						;
					} else {
						// x^y --> no go.
						throw unsuitable("Non-unity Trigonometric array is not suitable for positive real "
							"exponentiation.");
					}
				}
				return retval;
			}
			template <class ArgsTuple>
			Derived root(const max_fast_int &n, const ArgsTuple &args_tuple) const {
				if (n == 0) {
					throw division_by_zero();
				} else if (n == 1) {
					return Derived(*derived_const_cast);
				}
				return pow(1. / static_cast<double>(n), args_tuple);
			}
			// NOTE: here args_tuple must be the merge of the series undergoing the substitution and
			// the series used for the substitution.
			template <class RetSeries, class PosTuple, class SubCaches, class ArgsTuple>
			RetSeries sub(const PosTuple &pos_tuple, SubCaches &sub_caches, const ArgsTuple &args_tuple) const {
				typedef typename RetSeries::term_type ret_term_type;
				typedef typename ret_term_type::cf_type ret_cf_type;
				RetSeries retval;
				// If the argument is not present here, the return series will have one term consisting
				// of a unitary coefficient and this very trig_array.
				if (!pos_tuple.template get<Derived::position>().first) {
					retval = key_series_builder::template run<RetSeries>(*derived_const_cast, args_tuple);
				} else {
					const size_t pos = pos_tuple.template get<Derived::position>().second;
					const max_fast_int power = static_cast<max_fast_int>((*derived_const_cast)[pos]);
					p_assert(pos < derived_const_cast->size());
					Derived tmp_ta(*derived_const_cast);
					// Let's turn off the multiplier associated to the symbol we are substituting.
					tmp_ta[pos] = 0;
					// Build the orig_cos series.
					tmp_ta.flavour() = true;
					RetSeries orig_cos = key_series_builder::template run<RetSeries>(tmp_ta,args_tuple);
					// Build the orig_sin series.
					tmp_ta.flavour() = false;
					RetSeries orig_sin = key_series_builder::template run<RetSeries>(tmp_ta,args_tuple);
					p_assert(retval.empty());
					if (derived_const_cast->flavour()) {
						retval.add(orig_cos, args_tuple);
						retval.mult_by(
							sub_caches.template get<Derived::position>()[power].real(args_tuple),
						args_tuple);
						orig_sin.mult_by(
							sub_caches.template get<Derived::position>()[power].imag(args_tuple),
						args_tuple);
						retval.subtract(orig_sin, args_tuple);
					} else {
						retval.add(orig_sin, args_tuple);
						retval.mult_by(
							sub_caches.template get<Derived::position>()[power].real(args_tuple),
						args_tuple);
						orig_cos.mult_by(
							sub_caches.template get<Derived::position>()[power].imag(args_tuple),
						args_tuple);
						retval.add(orig_cos, args_tuple);
					}
				}
				return retval;
			}
			template <class RetSeries, class PosTuple, class SubCaches, class ArgsTuple>
			RetSeries ei_sub(const PosTuple &pos_tuple, SubCaches &sub_caches, const ArgsTuple &args_tuple) const {
				typedef typename RetSeries::term_type ret_term_type;
				typedef typename ret_term_type::cf_type ret_cf_type;
				RetSeries retval;
				// If the argument is not present here, the return series will have one term consisting
				// of a unitary coefficient and this very trig_array.
				if (!pos_tuple.template get<Derived::position>().first) {
					retval.insert(ret_term_type(ret_cf_type(static_cast<max_fast_int>(1), args_tuple),
						*derived_const_cast), args_tuple);
				} else {
					const size_t pos = pos_tuple.template get<Derived::position>().second;
					const max_fast_int power = static_cast<max_fast_int>((*derived_const_cast)[pos]);
					p_assert(pos < derived_const_cast->size());
					Derived tmp_ta(*derived_const_cast);
					RetSeries orig_cos, orig_sin;
					tmp_ta[pos] = 0;
					orig_cos.insert(ret_term_type(ret_cf_type(static_cast<max_fast_int>(1), args_tuple), tmp_ta),
						args_tuple);
					tmp_ta.flavour() = false;
					orig_sin.insert(ret_term_type(ret_cf_type(static_cast<max_fast_int>(1), args_tuple), tmp_ta),
						args_tuple);
					p_assert(retval.empty());
					if (derived_const_cast->flavour()) {
						retval.add(orig_cos, args_tuple);
						retval.mult_by(
							sub_caches.template get<Derived::position>()[power].real(args_tuple),
						args_tuple);
						orig_sin.mult_by(
							sub_caches.template get<Derived::position>()[power].imag(args_tuple),
						args_tuple);
						retval.subtract(orig_sin, args_tuple);
					} else {
						retval.add(orig_sin, args_tuple);
						retval.mult_by(
							sub_caches.template get<Derived::position>()[power].real(args_tuple),
						args_tuple);
						orig_cos.mult_by(
							sub_caches.template get<Derived::position>()[power].imag(args_tuple),
						args_tuple);
						retval.add(orig_cos, args_tuple);
					}
				}
				return retval;
			}
		protected:
			trig_array_commons() {}
			explicit trig_array_commons(const std::string &s) {
				typedef typename Derived::value_type value_type;
				std::vector<std::string> sd;
				boost::split(sd, s, boost::is_any_of(std::string(1, derived_const_cast->separator)));
				// TODO: check here that we are not loading too many multipliers, outside trig_size_t range.
				// TODO: do it everywhere!
				const size_t w = sd.size();
				if (w == 0) {
					// Set flavour to true, so that trig_array is logically equivalent to unity.
					derived_cast->m_flavour = true;
					return;
				}
				// Now we know  w >= 1.
				derived_cast->resize(w - 1);
				for (size_t i = 0;i < w - 1; ++i) {
					(*derived_cast)[i] = utils::lexical_converter<value_type>(sd[i]);
				}
				// Take care of flavour.
				if (*sd.back().c_str() == 's') {
					derived_cast->m_flavour = false;
				} else if (*sd.back().c_str() == 'c') {
					derived_cast->m_flavour = true;
				} else {
					std::cout << "Warning: undefined flavour '" << sd.back() <<
						"', defaulting to cosine." << std::endl;
					derived_cast->m_flavour = true;
				}
			}
		private:
			// NOTICE: is there some caching mechanism that can be used here?
			template <int N, class ArgsTuple>
			double combined_time_eval(const ArgsTuple &args_tuple) const {
				p_static_check(N >= 0, "");
				const size_t w = derived_const_cast->size();
				p_assert(w <= args_tuple.template get<Derived::position>().size());
				double retval = 0.;
				for (size_t i = 0;i < w;++i) {
					// We must be sure that there actually is component N in every symbol we are going to use.
					if (args_tuple.template get<Derived::position>()[i]->time_eval().size() > N) {
						retval += (*derived_const_cast)[i] *
							args_tuple.template get<Derived::position>()[i]->time_eval()[N];
					}
				}
				return retval;
			}
	};
}

#undef derived_const_cast
#undef derived_cast

#endif
