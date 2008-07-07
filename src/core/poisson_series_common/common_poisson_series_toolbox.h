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
			std::complex<Derived> complexp() const {
				typedef typename std::complex<Derived>::term_type complex_term_type;
				typedef typename Derived::const_sorted_iterator const_sorted_iterator;
				typedef typename Derived::term_type::cf_type::term_type::cf_type poly_cf_type;
				typedef typename std::complex<Derived>::term_type::cf_type complex_cf_type;
				// Get the term that has unity trig vector and whose coefficient is a linear polynomial with integer
				// coefficients or a linear polynomial with integer coefficients and a single coefficient.
				std::pair<const_sorted_iterator, std::pair<std::vector<poly_cf_type>, std::vector<max_fast_int> > >
				int_linear_term = get_int_linear_term<const_sorted_iterator, poly_cf_type>();
				// Expand using Jacobi-Anger's identity.
				std::complex<Derived> retval;
				jacang_ancestor::jacobi_anger(int_linear_term.first, retval, derived_const_cast->m_arguments);
				retval.m_arguments = derived_const_cast->m_arguments;
				// If the linear term was found, take care of it.
				if (int_linear_term.first != derived_const_cast->template nth_index<0>().end()) {
					std::complex<Derived> tmp_series;
					// Assign poly args of this as trig args of tmp_series.
					tmp_series.m_arguments.template get<1>() = derived_const_cast->m_arguments.template get<0>();
					// Let's build the term to be inserted in tmp_series.
					complex_term_type tmp_term1;
					complex_term_type tmp_term2;
					tmp_term2.m_key.flavour() = false;
					tmp_term1.m_cf = complex_cf_type(std::complex<max_fast_int>(1, 0), tmp_series.m_arguments);
					tmp_term2.m_cf = complex_cf_type(std::complex<max_fast_int>(0, 1), tmp_series.m_arguments);
					tmp_term1.m_key.assign_int_vector(int_linear_term.second.second);
					tmp_term2.m_key.assign_int_vector(int_linear_term.second.second);
					tmp_series.insert(tmp_term1, tmp_series.m_arguments);
					tmp_series.insert(tmp_term2, tmp_series.m_arguments);
					// Take care of the coefficient-only term, if any.
					if (int_linear_term.second.first.size() > 0) {
						std::complex<Derived> tmp_series2;
						complex_term_type tmp_term;
						tmp_term.m_cf.insert(
							typename std::complex<Derived>::term_type::cf_type::term_type(
								int_linear_term.second.first[0].complexp(tmp_series2.m_arguments),
								typename std::complex<Derived>::term_type::cf_type::term_type::key_type()
							),
							tmp_series2.m_arguments);
						tmp_series2.insert(tmp_term, tmp_series2.m_arguments);
						tmp_series *= tmp_series2;
					}
					retval *= tmp_series;
				}
				return retval;
			}
			Derived cos() const {
				return complexp().real();
			}
			Derived sin() const {
				return complexp().imag();
			}
		private:
			template <class Iterator, class PolyCf>
			std::pair<Iterator, std::pair<std::vector<PolyCf>, std::vector<max_fast_int> > > get_int_linear_term() const {
				BOOST_STATIC_ASSERT((boost::is_same<PolyCf, typename Derived::term_type::cf_type::term_type::cf_type>::value));
				typedef typename Derived::const_sorted_iterator const_sorted_iterator;
				BOOST_STATIC_ASSERT((boost::is_same<Iterator, const_sorted_iterator>::value));
				const const_sorted_iterator it_f = derived_const_cast->template nth_index<0>().end();
				std::pair<const_sorted_iterator, std::pair<std::vector<PolyCf>, std::vector<max_fast_int> > > retval;
				retval.first = it_f;
				// Make space to accomodate all the elements of the linear combination.
				retval.second.second.resize(derived_const_cast->m_arguments.template get<0>().size());
				for (const_sorted_iterator it = derived_const_cast->template nth_index<0>().begin(); it != it_f; ++it) {
					// If the term's trigonometric part is unity, let's see if we can extract a linear combination of arguments
					// from the corresponding polynomial.
					if (it->m_key.is_unity()) {
						try {
							it->m_cf.get_int_linear_combination(retval.second);
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
