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

#ifndef PIRANHA_BASE_SERIES_COMPLEX_TOOLBOX_H
#define PIRANHA_BASE_SERIES_COMPLEX_TOOLBOX_H

#include <complex>

#include "../exceptions.h"
#include "../p_assert.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	template <class RealDerived>
	class base_series_complex_toolbox
	{
			typedef std::complex<RealDerived> Derived;
		public:
			template <class ArgsTuple>
			RealDerived real() const
			{
				RealDerived retval;
				return retval;
			}
			template <class ArgsTuple>
			RealDerived imag() const
			{
				RealDerived retval;
				return retval;
			}
			template <class ArgsTuple>
			Derived &mult_by(const std::complex<max_fast_int> &cn, const ArgsTuple &args_tuple) {
				if (cn.real() == 0 and cn.imag() == 0) {
					Derived retval;
					derived_cast->swap_terms(retval);
				} else if (cn.real() != 1 or cn.imag() != 0) {
					Derived retval(derived_cast->multiply_coefficients_by(cn, args_tuple));
					derived_cast->swap_terms(retval);
				}
				return *derived_cast;
			}
			template <class ArgsTuple>
			Derived &mult_by(const std::complex<double> &cx, const ArgsTuple &args_tuple) {
				if (cx.real() == 0 and cx.imag() == 0) {
					Derived retval;
					derived_cast->swap_terms(retval);
				} else if (cx.real() = ! 1 or cx.imag() != 0) {
					Derived retval(derived_cast->multiply_coefficients_by(cx, args_tuple));
					derived_cast->swap_terms(retval);
				}
				return *derived_cast;
			}
			template <class ArgsTuple>
			Derived &mult_by(const RealDerived &r, const ArgsTuple &args_tuple)
			{
				Derived retval(derived_cast->multiply_by_series(r, args_tuple));
				// Grab the terms accumulated into return value.
				derived_cast->swap_terms(retval);
				return *derived_cast;
			}
			template <class ArgsTuple>
			Derived &divide_by(const std::complex<max_fast_int> &cn, const ArgsTuple &args_tuple) {
				if (cn.real() == 0 and cn.imag() == 0) {
					throw division_by_zero();
				} else if (cn.real() != 1 or cn.imag() != 0) {
					Derived retval(derived_cast->divide_coefficients_by(cn, args_tuple));
					derived_cast->swap_terms(retval);
				}
				return *derived_cast;
			}
			template <class ArgsTuple>
			Derived &divide_by(const std::complex<double> &cx, const ArgsTuple &args_tuple) {
				if (cx.real() == 0 and cx.imag() == 0) {
					throw division_by_zero();
				} else if (cx.real() != 1 or cx.imag() != 0) {
					Derived retval(derived_cast->divide_coefficients_by(cx, args_tuple));
					derived_cast->swap_terms(retval);
				}
				return *derived_cast;
			}
		protected:
			template <class ArgsTuple>
			void construct_from_real(const RealDerived &r, const ArgsTuple &args_tuple)
			{
				// TODO: maybe this can be named construct_from_other_series and placed into base_series?
				typedef typename RealDerived::const_sorted_iterator real_iterator;
				typedef typename Derived::const_sorted_iterator complex_iterator;
				// Make sure we are being called from an empty series.
				p_assert(derived_const_cast->template nth_index<0>().empty());
				const real_iterator r_it_f = r.template nth_index<0>().end();
				complex_iterator it_hint = derived_const_cast->template nth_index<0>().end();
				for (real_iterator r_it = r.template nth_index<0>().begin(); r_it != r_it_f; ++r_it)
				{
					it_hint = derived_cast->insert(*r_it,args_tuple,it_hint);
				}
			}
			template <class ArgsTuple>
			void construct_from_real_imag(const RealDerived &r, const RealDerived &i, const ArgsTuple &args_tuple)
			{
				typedef typename RealDerived::const_sorted_iterator real_iterator;
				typedef typename Derived::const_sorted_iterator complex_iterator;
				typedef typename Derived::term_type complex_term_type;
				// Make sure we are being called from an empty series.
				p_assert(derived_const_cast->template nth_index<0>().empty());
				// Let's build the real part first.
				construct_from_real(r,args_tuple);
				// Now let's proceed to the imaginary part.
				const real_iterator i_it_f = i.template nth_index<0>().end();
				complex_term_type tmp;
				complex_iterator it_hint = derived_const_cast->template nth_index<0>().end();
				for (real_iterator i_it = i.template nth_index<0>().begin(); i_it != i_it_f; ++i_it)
				{
					tmp.m_key = i_it->m_key;
					tmp.m_cf.imag(i_it->m_cf,args_tuple);
					it_hint = derived_cast->insert(tmp,args_tuple,it_hint);
				}
			}
	};
}

#undef derived_const_cast
#undef derived_cast

#endif
