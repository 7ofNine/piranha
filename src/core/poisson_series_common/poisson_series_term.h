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

#ifndef PIRANHA_POISSON_SERIES_TERM_H
#define PIRANHA_POISSON_SERIES_TERM_H

#include <boost/tuple/tuple.hpp>
#include <string>

#include "../base_classes/base_term.h"
#include "../ntuple.h"

namespace piranha
{
	/// Term class for Fourier series.
	template <class Cf, class Trig, char Separator, class Allocator>
	class poisson_series_term: public base_term<Cf, Trig, Separator, Allocator, poisson_series_term<Cf, Trig, Separator, Allocator> >
	{
			/// Alias for the ancestor.
			typedef base_term<Cf, Trig, Separator, Allocator, poisson_series_term<Cf, Trig, Separator, Allocator> > ancestor;
			/// Alias for evaluation type.
			typedef typename ancestor::eval_type eval_type;
		public:
			/// Alias for coefficient type.
			typedef Cf cf_type;
			/// Alias for trigonometric type.
			typedef Trig trig_type;
			/// Result of the multiplication of two terms.
			typedef typename boost::tuple<poisson_series_term, poisson_series_term> multiplication_result;
			/// Operators used during term multiplication.
			//typedef typename boost::tuples<subtraction,addition> multiplication_key_operators;
			/// Default constructor.
			explicit poisson_series_term(): ancestor() {}
			/// Ctor from string.
			template <class ArgsTuple>
			explicit poisson_series_term(const std::string &str, const ArgsTuple &args_tuple):
					ancestor(str, args_tuple) {}
			/// Constructor from coefficient and trigonometric part.
			explicit poisson_series_term(const cf_type &c, const trig_type &t): ancestor(c, t) {}
			/// Generic copy constructor.
			/**
			 * Constructs from piranha::poisson_series_term with optionally different coefficient type.
			 */
			template <class Cf2, class ArgsTuple>
			explicit poisson_series_term(const poisson_series_term<Cf2, Trig, Separator, Allocator> &term, const ArgsTuple &a):
					ancestor(term, a) {}
			/// Smarter numerical evaluation
			/**
			 * Similar to brute force evaluation, with the difference that sine and cosine of trigonometric arguments are cached
			 * and re-used over the evaluation of the series. Typically faster by a factor of 2-3, depending on
			 * the series' characteristics.
			 * @param[in] te piranha::trig_evaluator object that caches complex exponentials of trigonometric arguments.
			 */
			template <class TrigEvaluator>
			eval_type t_eval(TrigEvaluator &te) const {
				eval_type retval = ancestor::m_cf.t_eval(te.m_value, te.m_ps->m_arguments);
				retval *= ancestor::m_key.t_eval(te);
				return retval;
			}
			/// Check if the term is canonical.
			template <class ArgsTuple>
			bool is_canonical(const ArgsTuple &) const {
				return (ancestor::m_key.sign() > 0);
			}
			// TODO: check if it makes sense to skip the check here and assume canonicalise will be used iff
			// is_canonical has already been tested.
			/// Canonicalise the term.
			template <class ArgsTuple>
			void canonicalise(const ArgsTuple &args_tuple) {
				if (!is_canonical(args_tuple)) {
					invert_trig_args(args_tuple);
				}
			}
			/// Term multiplication.
			/**
			 * NOTE: the result of multiplication here _must_ be canonical.
			 */
			template <class Cf2, class ArgsTuple>
			static void multiply(const poisson_series_term *t1, const poisson_series_term<Cf2,trig_type,Separator,Allocator> *t2,
								 multiplication_result &res, const ArgsTuple &args_tuple) {
				// Perform the trigonometric multiplication.
				t1->m_key.multiply(t2->m_key, res.template get<0>().m_key, res.template get<1>().m_key);
				// Handle coefficient multiplication. Do the first coefficient, then assign the second one.
				// TODO: maybe provide the semantics to coefficients for something like this:
				// cf1.multiply_by_cf(cf2,res.template get<0>().m_cf,args_tuple),
				// so that we can avoid a copy.
				res.template get<0>().m_cf = t1->m_cf;
				res.template get<0>().m_cf.mult_by(t2->m_cf, args_tuple);
				res.template get<0>().m_cf.divide_by((max_fast_int)2, args_tuple);
				res.template get<1>().m_cf = res.template get<0>().m_cf;
				// Now adjust the signs according to werner's formulas.
				if (t1->m_key.flavour() == t2->m_key.flavour()) {
					res.template get<0>().m_key.flavour() = res.template get<1>().m_key.flavour() = true;
					if (!t1->m_key.flavour()) {
						res.template get<1>().m_cf.invert_sign(args_tuple);
					}
				} else {
					res.template get<0>().m_key.flavour() = res.template get<1>().m_key.flavour() = false;
					if (t1->m_key.flavour()) {
						res.template get<0>().m_cf.invert_sign(args_tuple);
					}
				}
				// Finally, canonicalise the retval terms.
				res.template get<0>().canonicalise(args_tuple);
				res.template get<1>().canonicalise(args_tuple);
			}
		private:
			// Invert the sign of trigonometric multipliers.
			template <class ArgsTuple>
			void invert_trig_args(const ArgsTuple &args_tuple) {
				ancestor::m_key.invert_sign();
				if (!(ancestor::m_key.flavour())) {
					ancestor::m_cf.invert_sign(args_tuple);
				}
			}
	};
}
#endif
