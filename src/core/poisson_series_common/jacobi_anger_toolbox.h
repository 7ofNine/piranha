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

#ifndef PIRANHA_JACOBI_ANGER_TOOLBOX_H
#define PIRANHA_JACOBI_ANGER_TOOLBOX_H

#include <boost/static_assert.hpp>
#include <boost/type_traits/is_same.hpp>
#include <complex>

#include "../integer_typedefs.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	template <int TrigPos, class Derived>
	class jacobi_anger_toolbox
	{
			BOOST_STATIC_ASSERT(TrigPos >= 0);
		protected:
			template <class Iterator, class ArgsTuple>
			std::complex<Derived> jacobi_anger(const Iterator &it_avoid, const ArgsTuple &args_tuple) const {
				typedef typename Derived::const_sorted_iterator const_sorted_iterator;
				BOOST_STATIC_ASSERT((boost::is_same<Iterator, const_sorted_iterator>::value));
				std::complex<Derived> retval((max_fast_int)1, args_tuple);
				if (derived_const_cast->empty()) {
					return retval;
				}
				// We want to proceed backwards here.
				const const_sorted_iterator it_f = derived_const_cast->template nth_index<0>().begin();
				const_sorted_iterator it = derived_const_cast->template nth_index<0>().end();
				while (it != it_f) {
					--it;
					// Skip the iterator we want to avoid.
					if (it != it_avoid) {
						retval.mult_by(jacang_term(it, args_tuple), args_tuple);
					}
				}
				return retval;
			}
		private:
			template <class Iterator, class ArgsTuple>
			static std::complex<Derived> jacang_term(const Iterator &it, const ArgsTuple &args_tuple) {
				typedef typename std::complex<Derived>::term_type complex_term_type;
				// Let's determine the limit of the Jacobi-Anger development from the truncator of the series.
				// The Jacobi-Anger development is a development into bessel functions of the first kind starting
				// from zero and increasing in unity steps, hence the power_series_limit function can be used
				// straightforwardly.
				const size_t n = Derived::multiplier_type::truncator_type::power_series_limit(it->m_cf, args_tuple);
				std::complex<Derived> retval;
				{
					complex_term_type tmp_term;
					tmp_term.m_cf.real(it->m_cf.besselJ(0, args_tuple), args_tuple);
					retval.insert(tmp_term, args_tuple);
				}
				const size_t w = args_tuple.template get<TrigPos>().size();
				std::vector<max_fast_int> tmp_trig_mults(w);
				std::complex<max_fast_int> cos_multiplier(0, 2);
				for (size_t i = 1; i <= n; ++i) {
					complex_term_type tmp_term;
					tmp_term.m_cf.real(it->m_cf.besselJ((max_fast_int)i, args_tuple), args_tuple);
					it->m_key.upload_ints_to(tmp_trig_mults);
					for (size_t j = 0; j < w; ++j) {
						tmp_trig_mults[j] *= i;
					}
					tmp_term.m_key.assign_int_vector(tmp_trig_mults);
					switch (it->m_key.flavour()) {
					case true:
						tmp_term.m_cf.mult_by(cos_multiplier, args_tuple);
						break;
					case false:
						if (i % 2 == 0) {
							tmp_term.m_cf.mult_by((max_fast_int)2, args_tuple);
						} else {
							tmp_term.m_cf.mult_by(std::complex<max_fast_int>(0, 2), args_tuple);
							tmp_term.m_key.flavour() = false;
						}
					}
					retval.insert(tmp_term, args_tuple);
					// Update the multiplier for cosine terms.
					cos_multiplier *= std::complex<max_fast_int>(0, 1);
				}
				return retval;
			}

	};
}

#undef derived_const_cast
#undef derived_cast

#endif
