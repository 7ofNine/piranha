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

#ifndef PIRANHA_COMMON_POISSON_SERIES_TOOLBOX_H
#define PIRANHA_COMMON_POISSON_SERIES_TOOLBOX_H

#include <boost/static_assert.hpp>
#include <boost/type_traits/is_same.hpp>
#include <complex>
#include <utility>
#include <vector>

#include "../exceptions.h"
#include "../integer_typedefs.h"
#include "../p_assert.h"
#include "jacobi_anger_toolbox.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	template <class Derived>
	class common_poisson_series_toolbox: public jacobi_anger_toolbox<1, Derived>
	{
			typedef jacobi_anger_toolbox<1, Derived> jacang_ancestor;
		public:
			// NOTICE: this method assumes that the input args tuple already hase merged in as
			// trig arguments the poly arguments (see also below).
			template <class ArgsTuple>
			std::complex<Derived> complexp(const ArgsTuple &args_tuple) const {
				typedef typename std::complex<Derived>::term_type complex_term_type;
				typedef typename Derived::const_sorted_iterator const_sorted_iterator;
				typedef typename Derived::term_type::cf_type::term_type::cf_type poly_cf_type;
				typedef typename std::complex<Derived>::term_type::cf_type complex_cf_type;
				// Get the term that has unity trig vector and whose coefficient is a linear polynomial with integer
				// coefficients or a linear polynomial with integer coefficients and a single coefficient.
				std::pair<const_sorted_iterator, std::pair<std::vector<poly_cf_type>, std::vector<max_fast_int> > >
				int_linear_term(get_int_linear_term<const_sorted_iterator, poly_cf_type>(args_tuple));
				// Expand using Jacobi-Anger's identity.
				std::complex<Derived> retval;
				jacang_ancestor::jacobi_anger(int_linear_term.first, retval, args_tuple);
				// If the linear term was found, take care of it.
				if (int_linear_term.first != derived_const_cast->template nth_index<0>().end()) {
					std::complex<Derived> tmp_series;
					// Let's build the term to be inserted in tmp_series.
					complex_term_type tmp_term1;
					complex_term_type tmp_term2;
					tmp_term1.m_cf = complex_cf_type(std::complex<max_fast_int>(1, 0), args_tuple);
					tmp_term1.m_key.assign_int_vector(int_linear_term.second.second);
					tmp_term1.m_key.flavour() = true;
					tmp_term2.m_cf = complex_cf_type(std::complex<max_fast_int>(0, 1), args_tuple);
					tmp_term2.m_key.assign_int_vector(int_linear_term.second.second);
					tmp_term2.m_key.flavour() = false;
					tmp_series.insert(tmp_term1, args_tuple);
					tmp_series.insert(tmp_term2, args_tuple);
					// Take care of the numerical-coefficient-only term, if any.
					if (int_linear_term.second.first.size() > 0) {
						std::complex<Derived> tmp_series2;
						// NOTE: tmp_series2's arguments tuple is empty, it does not matter
						// because we are explicitly dealing with single numerical coefficient series.
						complex_term_type tmp_term;
						tmp_term.m_cf.insert(
							typename std::complex<Derived>::term_type::cf_type::term_type(
								int_linear_term.second.first[0].complexp(tmp_series2.m_arguments),
								typename std::complex<Derived>::term_type::cf_type::term_type::key_type()
							),
							tmp_series2.m_arguments);
						tmp_series2.insert(tmp_term, tmp_series2.m_arguments);
						tmp_series.mult_by(tmp_series2,args_tuple);
					}
					retval.mult_by(tmp_series,args_tuple);
				}
				return retval;
			}
			std::complex<Derived> complexp() const {
				// In order to account for a potential integer linear combination of arguments
				// we must merge in as trigonometric arguments the polynomial arguments. The safe
				// way to do this is using named_series::merge_args with a phony series having zero
				// polynomial arguments and as trigonometric arguments the polynomial arguments of this.
				Derived copy(*derived_const_cast), tmp;
				tmp.m_arguments.template get<1>() = derived_const_cast->m_arguments.template get<0>();
				copy.merge_args(tmp);
				// Now we can call in the complexp method from above.
				std::complex<Derived> retval(copy.complexp(copy.m_arguments));
				retval.m_arguments = copy.m_arguments;
				retval.trim();
				return retval;
			}
			Derived cos() const {
				return complexp().real();
			}
			Derived sin() const {
				return complexp().imag();
			}
		private:
			template <class Iterator, class PolyCf, class ArgsTuple>
			std::pair<Iterator, std::pair<std::vector<PolyCf>, std::vector<max_fast_int> > >
			get_int_linear_term(const ArgsTuple &args_tuple) const {
				BOOST_STATIC_ASSERT((boost::is_same<PolyCf, typename Derived::term_type::cf_type::term_type::cf_type>::value));
				typedef typename Derived::const_sorted_iterator const_sorted_iterator;
				BOOST_STATIC_ASSERT((boost::is_same<Iterator, const_sorted_iterator>::value));
				const const_sorted_iterator it_f = derived_const_cast->template nth_index<0>().end();
				std::pair<const_sorted_iterator, std::pair<std::vector<PolyCf>, std::vector<max_fast_int> > > retval;
				retval.first = it_f;
				// Make space to accomodate all the elements of the linear combination.
				// We need as much space as the number ofr trig args.
				retval.second.second.resize(args_tuple.template get<1>().size());
				for (const_sorted_iterator it = derived_const_cast->template nth_index<0>().begin(); it != it_f; ++it) {
					// If the term's trigonometric part is unity, let's see if we can extract a linear combination of arguments
					// from the corresponding polynomial.
					if (it->m_key.is_unity()) {
						try {
							it->m_cf.template get_int_linear_combination<1>(retval.second,args_tuple);
						} catch (const unsuitable &) {
							// If we are unable to extract a proper linear combination from the unity term, erase retval
							// and break out.
							retval.second.first.clear();
							retval.second.second.clear();
							break;
						}
						retval.first = it;
						break;
					}
				}
				p_assert(retval.second.first.size() <= 1);
				return retval;
			}
	};
}

#undef derived_const_cast
#undef derived_cast

#endif
