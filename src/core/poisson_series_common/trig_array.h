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

#ifndef PIRANHA_TRIG_ARRAY_H
#define PIRANHA_TRIG_ARRAY_H

#include <boost/algorithm/string/split.hpp>
#include <cmath> // For std::abs.
#include <complex> // For std::complex<SubSeries>.
#include <iostream>
#include <memory> // For standard allocator.
#include <string>
#include <utility> // For std::pair.
#include <vector>

#include "../base_classes/int_array.h"
#include "../base_classes/series_builders.h"
#include "../base_classes/toolbox.h"
#include "../common_functors.h"
#include "../config.h"
#include "../exceptions.h"
#include "../int_power_cache.h"
#include "../psym.h"
#include "../settings.h"
#include "../utils.h" // For lexical converter.
// TODO: remove when efficient trig evaluation is implemented.
#include "trig_evaluator.h"

#define __PIRANHA_TRIG_ARRAY_TP_DECL int Bits, int Pos, class Allocator
#define __PIRANHA_TRIG_ARRAY_TP Bits,Pos,Allocator

namespace piranha
{
	template < __PIRANHA_TRIG_ARRAY_TP_DECL = std::allocator<char> >
	struct trig_array {
		typedef toolbox<trig_array<__PIRANHA_TRIG_ARRAY_TP> > type;
	};

	/// Trigonometric array, dynamically sized version.
	/**
	 * It wraps a piranha::int_array with signed integer sized Bits, and adds the
	 * capabilities needed for trigonometric manipulation.
	 */
	template < __PIRANHA_TRIG_ARRAY_TP_DECL>
	class toolbox<trig_array<__PIRANHA_TRIG_ARRAY_TP> >: public int_array<__PIRANHA_TRIG_ARRAY_TP, toolbox<trig_array<__PIRANHA_TRIG_ARRAY_TP> > >
	{
			typedef int_array<__PIRANHA_TRIG_ARRAY_TP, toolbox<trig_array<__PIRANHA_TRIG_ARRAY_TP> > > ancestor;
			friend class int_array<__PIRANHA_TRIG_ARRAY_TP, toolbox<trig_array<__PIRANHA_TRIG_ARRAY_TP> > >;
			template <class SubSeries, class ArgsTuple>
			class sub_cache: public int_power_cache<std::complex<SubSeries>,
				typename base_series_arithmetics<std::complex<SubSeries>,ArgsTuple>::type>
			{
					typedef int_power_cache<std::complex<SubSeries>,
						typename base_series_arithmetics<std::complex<SubSeries>,ArgsTuple>::type> ancestor;
					enum status {
						zero,
						one,
						full
					};
				public:
					sub_cache():ancestor::int_power_cache(),m_status(zero),m_errmsg() {}
					void setup(const SubSeries &s, const ArgsTuple *args_tuple) {
						this->m_arith_functor.m_args_tuple = args_tuple;
						this->m_container[0] = std::complex<SubSeries>().base_add(1,*args_tuple);
						try {
							std::complex<SubSeries> tmp1 = s.ei(*args_tuple);
							this->m_container[1] = tmp1;
							m_status = one;
							SubSeries tmp2(s);
							tmp2.base_mult_by(-1,*args_tuple);
							std::complex<SubSeries> tmp3 = tmp2.ei(*args_tuple);
							this->m_container[-1] = tmp3;
							m_status = full;
						} catch (const unsuitable &u) {
							m_errmsg = u.what();
						}
					}
					const std::complex<SubSeries> &operator[](const int &n) {
						switch (m_status) {
							case zero:
								if (n != 0) {
									throw unsuitable(std::string("The substitution cache was unable to "
										"compute the complex exponential of the series used for substitution. "
										"The reported error was:\n") + m_errmsg);
								}
							case one:
								if (n < 0) {
									throw unsuitable(std::string("The substitution cache was unable to "
										"compute the inverse complex exponential of the series used for substitution. "
										"The reported error was:\n") + m_errmsg);
								}
							default:
								;
						}
						return ancestor::operator[](n);
					}
				private:
					status		m_status;
					std::string	m_errmsg;
			};
			template <class SubSeries, class ArgsTuple>
			class ei_sub_cache: public int_power_cache<SubSeries, typename base_series_arithmetics<SubSeries,ArgsTuple>::type>
			{
					typedef int_power_cache<SubSeries, typename base_series_arithmetics<SubSeries,ArgsTuple>::type> ancestor;
				public:
					ei_sub_cache():ancestor::int_power_cache() {}
					void setup(const SubSeries &s, const ArgsTuple *args_tuple) {
						this->m_arith_functor.m_args_tuple = args_tuple;
						this->m_container[0] = SubSeries().base_add(1,*args_tuple);
						this->m_container[1] = s;
					}
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
				// TODO: check here that we are not loading too many multipliers, outside trig_size_t range.
				// TODO: do it everywhere!
				const size_t w = sd.size();
				if (w == 0) {
					// Set flavour to true, so that trig_array is logically equivalent to unity.
					this->m_flavour = true;
					return;
				}
				// Now we know  w >= 1.
				this->resize(w - 1);
				for (size_t i = 0;i < w - 1; ++i) {
					(*this)[i] = utils::lexical_converter<value_type>(sd[i]);
				}
				// Take care of flavour.
				if (*sd.back().c_str() == 's') {
					this->m_flavour = false;
				} else if (*sd.back().c_str() == 'c') {
					this->m_flavour = true;
				} else {
					std::cout << "Warning: undefined flavour '" << sd.back() <<
						"', defaulting to cosine." << std::endl;
					this->m_flavour = true;
				}
			}
			template <class ArgsTuple>
			explicit toolbox(const psym &p, const int &n, const ArgsTuple &a): ancestor::int_array(p, n, a) {}
			template <int Pos2>
			explicit toolbox(const toolbox<trig_array<Bits,Pos2,Allocator> > &ta): ancestor::int_array(ta) {}
			// Probing.
			/// Data footprint.
			/**
			 * Returns the memory occupied by the data members.
			 */
			size_t data_footprint() const {
				return (ancestor::size()*sizeof(value_type));
			}
			// Math.
			/// Multiplication.
			/**
			 * Used in poisson_series_term multiplication.
			 * TODO: update docs below.
			 * Multiplication of two trigonometric functions using Werner's formulas, i.e.
			 * \f[
			 * C\cos\alpha\cdot\cos\beta=
			 * \frac{C}{2} \cos \left( \alpha - \beta \right) + \frac{C}{2} \cos \left( \alpha + \beta \right)
			 * \f]
			 * and the likes. Notice that in the first return value always goes the \f$ \alpha - \beta \f$ term
			 * and in the second one always goes \f$ \alpha + \beta \f$ one.
			 * Please also note that no assumptions are made with respect to return values' content
			 * (e.g., it is not guaranteed that return values are empty).
			 * @param[in] t2 factor.
			 * @param[out] ret1 first return value.
			 * @param[out] ret2 second return value.
			 */
			template <class TrigArray>
			void multiply(const TrigArray &t2, toolbox &ret1, toolbox &ret2) const
			// NOTE: we are not using here a general version of vector addition/subtraction
			// because this way we can do two operations (+ and -) every cycle. This is a performance
			// critical part, so the optimization should be worth the hassle.
			{
				const size_type max_w = ancestor::size(), min_w = t2.size();
				// Assert widths, *this should always come from a regular Poisson series, and its width should hence be
				// already adjusted my merge_args in multiplication routines.
				p_assert(max_w >= min_w);
				// Adjust the width of retvals, if needed.
				ret1.resize(max_w);
				ret2.resize(max_w);
				p_assert(ret1.size() == max_w);
				p_assert(ret2.size() == max_w);
				size_type i;
				for (i = 0; i < min_w; ++i) {
					ret1[i] = (*this)[i] - t2[i];
					ret2[i] = (*this)[i] + t2[i];
				}
				for (; i < max_w; ++i) {
					ret1[i] = (*this)[i];
					ret2[i] = (*this)[i];
				}
			}
			/// Get flavour.
			bool get_flavour() const {
				return this->m_flavour;
			}
			/// Set flavour.
			void set_flavour(bool f) {
				this->m_flavour = f;
			}
			// I/O.
			template <class ArgsTuple>
			void print_plain(std::ostream &out_stream, const ArgsTuple &args_tuple) const {
				// We assert like this because we want to make sure we don't go out of boundaries,
				// and because in case of fixed-width we may have smaller size of v wrt to "real" size.
				p_assert(args_tuple.template get<ancestor::position>().size() <= this->size());
				(void)args_tuple;
				this->print_elements(out_stream);
				// Print the separator before flavour only if we actually printed something above.
				if (this->size() != 0) {
					out_stream << this->separator;
				}
				if (this->m_flavour) {
					out_stream << "c";
				} else {
					out_stream << "s";
				}
			}
			template <class ArgsTuple>
			void print_pretty(std::ostream &out_stream, const ArgsTuple &args_tuple) const {
				p_assert(args_tuple.template get<ancestor::position>().size() <= this->m_size);
				if (this->m_flavour) {
					out_stream << "cos(";
				} else {
					out_stream << "sin(";
				}
				bool printed_something = false;
				for (size_t i = 0; i < this->m_size; ++i) {
					const int n = this->m_container.v[i];
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
						out_stream << args_tuple.template get<ancestor::position>()[i].get_name();
						printed_something = true;
					}
				}
				out_stream << ')';
			}
			void print_latex(std::ostream &out_stream, const vector_psym &v) const {
				const size_t w = v.size();
				p_assert(w <= this->size())
				switch (this->m_flavour) {
					case true:
						out_stream << "c&";
						break;
					case false:
						out_stream << "s&";
				}
				bool first_one = true;
				std::string tmp("$");
				for (size_t i = 0;i < w;++i) {
					if ((*this)[i] != 0) {
						if ((*this)[i] > 0 && !first_one) {
							tmp.append("+");
						}
						if ((*this)[i] == -1) {
							tmp.append("-");
						} else if ((*this)[i] == 1) {} else {
							tmp.append(boost::lexical_cast<std::string>((int)(*this)[i]));
						}
						tmp.append(v[i].get_name());
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
				return (this->elements_are_zero() && this->m_flavour);
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
				p_assert(args_tuple.template get<ancestor::position>().size() >= this->size());
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
				const size_t w = this->size();
				p_assert(w <= args_tuple.template get<ancestor::position>().size());
				double retval = 0.;
				for (size_t i = 0;i < w;++i) {
					if ((*this)[i] != 0) {
						retval += (*this)[i] * args_tuple.template get<ancestor::position>()[i].eval(t);
					}
				}
				if (this->m_flavour) {
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
				p_assert(w <= this->size());
				std::complex<double> retval(1.);
				for (size_t i = 0;i < w;++i) {
					if ((*this)[i] != 0) {
						retval *= te.request_ei(i, (*this)[i]);
					}
				}
				if (this->m_flavour) {
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
				const size_t w = this->size();
				for (size_t i = 0;i < w;++i) {
					if ((*this)[i] > 0) {
						return 1;
					}
					if ((*this)[i] < 0) {
						return -1;
					}
				}
				return 1;
			}
			// All multipliers are zero and flavour is sine.
			template <class ArgsTuple>
			bool is_ignorable(const ArgsTuple &) const {
				return (this->elements_are_zero() && !this->m_flavour);
			}
			/// Equality test.
			bool operator==(const toolbox &t2) const {
				return (this->m_flavour == t2.m_flavour && this->elements_equal_to(t2));
			}
			/// Less than.
			bool operator<(const toolbox &t2) const {
				if (this->m_flavour < t2.m_flavour) {
					return true;
				} else if (this->m_flavour > t2.m_flavour) {
					return false;
				}
				return this->lex_comparison(t2);
			}
			/// Calculate hash_value.
			/**
			 * Used by the hash_value overload for piranha::base_term.
			 */
			size_t hash_value() const {
				size_t retval = this->elements_hasher();
				boost::hash_combine(retval, this->m_flavour);
				return retval;
			}
			/// Multiply by an integer.
			void mult_by_int(const int &n) {
				const size_t w = this->size();
				for (size_t i = 0;i < w;++i) {
					(*this)[i] *= n;
				}
			}
			/// Calculate partial derivative.
			/**
			 * Result is a pair consisting of an integer and a trigonometric array.
			 */
			template <class Series, class PosTuple, class ArgsTuple>
			Series partial(const PosTuple &pos_tuple, const ArgsTuple &args_tuple) const {
				Series retval;
				// Do something only if the argument of the partial derivation is present in the trigonometric array.
				// Otherwise the above retval will return, and it will deliver a zero integer multiplier to be
				// multiplied by the coefficient in the partial derivation of the whole term.
				p_assert(pos_tuple.template get<ancestor::position>().size() == 1);
				if (pos_tuple.template get<ancestor::position>()[0].first) {
					toolbox copy(*this);
					const size_t pos = pos_tuple.template get<ancestor::position>()[0].second;
					// Change the flavour of the resulting key.
					copy.m_flavour = !this->m_flavour;
					p_assert(pos < this->size());
					retval = Series::base_series_from_key(copy,args_tuple);
					if (this->m_flavour) {
						retval.base_mult_by(-1,args_tuple);
					}
					retval.base_mult_by((*this)[pos],args_tuple);
				}
				return retval;
			}
			/// Exponentiation.
			template <class ArgsTuple>
			toolbox pow(const double &y, const ArgsTuple &) const {
				return pow_number(y);
			}
			template <class ArgsTuple>
			toolbox root(const int &n, const ArgsTuple &args_tuple) const {
				if (n == 0) {
					throw division_by_zero();
				} else if (n == 1) {
					return toolbox(*this);
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
				// NOTE: for now we can substitute one symbol at a time.
				p_assert(pos_tuple.template get<ancestor::position>().size() == 1);
				if (!pos_tuple.template get<ancestor::position>()[0].first) {
					retval = key_series_builder::template run<RetSeries>(*this, args_tuple);
				} else {
					const size_t pos = pos_tuple.template get<ancestor::position>()[0].second;
					const int power = static_cast<int>((*this)[pos]);
					p_assert(pos < this->size());
					toolbox tmp_ta(*this);
					// Let's turn off the multiplier associated to the symbol we are substituting.
					tmp_ta[pos] = 0;
					// NOTE: important: we need key builders here because we may be building RetSeries
					// whose key is _not_ a trig_array, in principle, so we cannot build a term consisting
					// of a trig_array and unity coefficient and simply insert it.
					// Build the orig_cos series.
					tmp_ta.set_flavour(true);
					RetSeries orig_cos = key_series_builder::template run<RetSeries>(tmp_ta,args_tuple);
					// Build the orig_sin series.
					tmp_ta.set_flavour(false);
					RetSeries orig_sin = key_series_builder::template run<RetSeries>(tmp_ta,args_tuple);
					p_assert(retval.empty());
					if (this->get_flavour()) {
						retval.base_add(orig_cos, args_tuple);
						retval.base_mult_by(
							sub_caches.template get<ancestor::position>()[power].base_real(args_tuple),
						args_tuple);
						orig_sin.base_mult_by(
							sub_caches.template get<ancestor::position>()[power].base_imag(args_tuple),
						args_tuple);
						retval.base_subtract(orig_sin, args_tuple);
					} else {
						retval.base_add(orig_sin, args_tuple);
						retval.base_mult_by(
							sub_caches.template get<ancestor::position>()[power].base_real(args_tuple),
						args_tuple);
						orig_cos.base_mult_by(
							sub_caches.template get<ancestor::position>()[power].base_imag(args_tuple),
						args_tuple);
						// NOTE: series multadd here (and multiply by -1 to do subtraction too)?
						// Below too...
						retval.base_add(orig_cos, args_tuple);
					}
				}
				return retval;
			}
			template <class RetSeries, class PosTuple, class SubCaches, class ArgsTuple>
			RetSeries ei_sub(const PosTuple &pos_tuple, SubCaches &sub_caches, const ArgsTuple &args_tuple) const {
				typedef typename RetSeries::term_type ret_term_type;
				typedef typename ret_term_type::cf_type ret_cf_type;
				RetSeries retval;
				p_assert(pos_tuple.template get<ancestor::position>().size() == 1);
				if (!pos_tuple.template get<ancestor::position>()[0].first) {
					retval = key_series_builder::template run<RetSeries>(*this, args_tuple);
				} else {
					const size_t pos = pos_tuple.template get<ancestor::position>()[0].second;
					const int power = static_cast<int>((*this)[pos]);
					p_assert(pos < this->size());
					toolbox tmp_ta(*this);
					tmp_ta[pos] = 0;
					tmp_ta.set_flavour(true);
					RetSeries orig_cos = key_series_builder::template run<RetSeries>(tmp_ta,args_tuple);
					tmp_ta.set_flavour(false);
					RetSeries orig_sin = key_series_builder::template run<RetSeries>(tmp_ta,args_tuple);
					p_assert(retval.empty());
					if (this->m_flavour) {
						retval.base_add(orig_cos, args_tuple);
						retval.base_mult_by(
							sub_caches.template get<ancestor::position>()[power].base_real(args_tuple),
						args_tuple);
						orig_sin.base_mult_by(
							sub_caches.template get<ancestor::position>()[power].base_imag(args_tuple),
						args_tuple);
						retval.base_subtract(orig_sin, args_tuple);
					} else {
						retval.base_add(orig_sin, args_tuple);
						retval.base_mult_by(
							sub_caches.template get<ancestor::position>()[power].base_real(args_tuple),
						args_tuple);
						orig_cos.base_mult_by(
							sub_caches.template get<ancestor::position>()[power].base_imag(args_tuple),
						args_tuple);
						retval.base_add(orig_cos, args_tuple);
					}
				}
				return retval;
			}
		private:
			// NOTICE: is there some caching mechanism that can be used here?
			template <int N, class ArgsTuple>
			double combined_time_eval(const ArgsTuple &args_tuple) const {
				p_static_check(N >= 0, "");
				const size_t w = this->size();
				p_assert(w <= args_tuple.template get<ancestor::position>().size());
				double retval = 0.;
				for (size_t i = 0;i < w;++i) {
					// We must be sure that there actually is component N in every symbol we are going to use.
					if (args_tuple.template get<ancestor::position>()[i].get_time_eval().size() > N) {
						retval += (*this)[i] *
							args_tuple.template get<ancestor::position>()[i].get_time_eval()[N];
					}
				}
				return retval;
			}
			template <class Number>
			toolbox pow_number(const Number &y) const {
				const bool int_zero = this->elements_are_zero();
				toolbox retval;
				if (y < 0) {
					if (int_zero && !this->m_flavour) {
						// 0**-y.
						throw division_by_zero();
					} else if (int_zero && this->m_flavour) {
						// 1**-y == 1. Don't do anything because retval is already initialized properly.
						;
					} else {
						// x**-y -> no go.
						throw unsuitable("Non-unity Trigonometric array is not suitable for negative exponentiation.");
					}
				} else if (y == 0) {
					// x**0 == 1. Don't do nothing because retval is already initialized properly.
					;
				} else {
					if (int_zero && !this->m_flavour) {
						// 0**y == 0.
						retval.m_flavour = false;
					} else if (int_zero && this->m_flavour) {
						// 1**y == 1. Don't do anything because retval is already initialized properly.
						;
					} else {
						// x**y --> no go.
						throw unsuitable("Non-unity Trigonometric array is not suitable for positive exponentiation.");
					}
				}
				return retval;
			}
	};
}

#undef __PIRANHA_TRIG_ARRAY_TP_DECL
#undef __PIRANHA_TRIG_ARRAY_TP

#endif
