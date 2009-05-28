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

#ifndef PIRANHA_NAMED_SERIES_MATH_H
#define PIRANHA_NAMED_SERIES_MATH_H

#include <vector>

#include "../exceptions.h"
#include "../mp.h"
#include "../ntuple.h"
#include "../psym.h"
#include "../settings.h"
#include "named_series_def.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	template <bool Sign, class Derived2>
	inline Derived &toolbox<named_series<__PIRANHA_NAMED_SERIES_TP> >::merge_with_series(const Derived2 &s2)
	{
		// If we are merging with self, create a copy and call recursively.
		if ((void *)derived_cast == (void *)(&s2)) {
			__PDEBUG(std::cout << "Merging with self, performing a copy." << '\n');
			merge_with_series<Sign>(Derived(*derived_const_cast));
		} else {
			merge_args(s2);
			derived_cast->template merge_terms<Sign>(s2, m_arguments);
		}
		trim();
		return *derived_cast;
	}

	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	template <class Derived2>
	inline Derived &toolbox<named_series<__PIRANHA_NAMED_SERIES_TP> >::mult_by_series(const Derived2 &s2)
	{
		// First we merge the arguments of the two series.
		merge_args(s2);
		// Then we perform the multiplication.
		derived_cast->base_mult_by(s2, m_arguments);
		trim();
		return *derived_cast;
	}

	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	template <bool Sign, class Number>
	inline Derived &toolbox<named_series<__PIRANHA_NAMED_SERIES_TP> >::merge_number_helper(const Number &x)
	{
		derived_cast->template merge_with_number<Sign>(x, m_arguments);
		trim();
		return *derived_cast;
	}

	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline Derived &toolbox<named_series<__PIRANHA_NAMED_SERIES_TP> >::operator+=(const double &x)
	{
		return merge_number_helper<true>(x);
	}

	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline Derived &toolbox<named_series<__PIRANHA_NAMED_SERIES_TP> >::operator+=(const mp_rational &q)
	{
		return merge_number_helper<true>(q);
	}

	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline Derived &toolbox<named_series<__PIRANHA_NAMED_SERIES_TP> >::operator+=(const mp_integer &z)
	{
		return merge_number_helper<true>(z);
	}

	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline Derived &toolbox<named_series<__PIRANHA_NAMED_SERIES_TP> >::operator+=(const Derived &s2)
	{
		return merge_with_series<true>(s2);
	}

	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline Derived &toolbox<named_series<__PIRANHA_NAMED_SERIES_TP> >::operator-=(const double &x)
	{
		return merge_number_helper<false>(x);
	}

	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline Derived &toolbox<named_series<__PIRANHA_NAMED_SERIES_TP> >::operator-=(const mp_rational &q)
	{
		return merge_number_helper<false>(q);
	}

	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline Derived &toolbox<named_series<__PIRANHA_NAMED_SERIES_TP> >::operator-=(const mp_integer &z)
	{
		return merge_number_helper<false>(z);
	}

	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline Derived &toolbox<named_series<__PIRANHA_NAMED_SERIES_TP> >::operator-=(const Derived &s2)
	{
		return merge_with_series<false>(s2);
	}

	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline Derived toolbox<named_series<__PIRANHA_NAMED_SERIES_TP> >::operator-() const
	{
		Derived retval(*derived_const_cast);
		retval *= -1;
		return retval;
	}

	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	template <class Number>
	inline Derived &toolbox<named_series<__PIRANHA_NAMED_SERIES_TP> >::mult_number_helper(const Number &x)
	{
		derived_cast->base_mult_by(x, m_arguments);
		trim();
		return *derived_cast;
	}

	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline Derived &toolbox<named_series<__PIRANHA_NAMED_SERIES_TP> >::operator*=(const double &x)
	{
		return mult_number_helper(x);
	}

	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline Derived &toolbox<named_series<__PIRANHA_NAMED_SERIES_TP> >::operator*=(const mp_rational &q)
	{
		return mult_number_helper(q);
	}

	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline Derived &toolbox<named_series<__PIRANHA_NAMED_SERIES_TP> >::operator*=(const mp_integer &z)
	{
		return mult_number_helper(z);
	}

	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline Derived &toolbox<named_series<__PIRANHA_NAMED_SERIES_TP> >::operator*=(const Derived &s2)
	{
		return mult_by_series(s2);
	}

	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	template <class Number>
	inline Derived &toolbox<named_series<__PIRANHA_NAMED_SERIES_TP> >::divide_number_helper(const Number &x)
	{
		derived_cast->base_divide_by(x, m_arguments);
		trim();
		return *derived_cast;
	}

	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline Derived &toolbox<named_series<__PIRANHA_NAMED_SERIES_TP> >::operator/=(const double &x)
	{
		return divide_number_helper(x);
	}

	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline Derived &toolbox<named_series<__PIRANHA_NAMED_SERIES_TP> >::operator/=(const mp_rational &q)
	{
		return divide_number_helper(q);
	}

	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline Derived &toolbox<named_series<__PIRANHA_NAMED_SERIES_TP> >::operator/=(const mp_integer &z)
	{
		return divide_number_helper(z);
	}

	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline Derived toolbox<named_series<__PIRANHA_NAMED_SERIES_TP> >::pow(const double &x) const
	{
		Derived retval(derived_const_cast->base_pow(x, m_arguments));
		retval.m_arguments = m_arguments;
		retval.trim();
		return retval;
	}

	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline Derived toolbox<named_series<__PIRANHA_NAMED_SERIES_TP> >::pow(const mp_rational &q) const
	{
		Derived retval(derived_const_cast->base_pow(q, m_arguments));
		retval.m_arguments = m_arguments;
		retval.trim();
		return retval;
	}

	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline Derived toolbox<named_series<__PIRANHA_NAMED_SERIES_TP> >::root(const int &n) const
	{
		Derived retval(derived_const_cast->base_root(n, m_arguments));
		retval.m_arguments = m_arguments;
		retval.trim();
		return retval;
	}

	/// Partial derivative with respect to a piranha::psym.
	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline Derived toolbox<named_series<__PIRANHA_NAMED_SERIES_TP> >::partial(const psym &p, const int &n) const
	{
		typedef typename ntuple<std::vector<std::pair<bool, size_t> >, n_arguments_sets>::type pos_tuple_type;
		const pos_tuple_type pos_tuple = psyms2pos(vector_psym(1,p),m_arguments);
		Derived retval(derived_const_cast->base_partial(n, pos_tuple, m_arguments));
		retval.m_arguments = m_arguments;
		retval.trim();
		return retval;
	}
}

#undef derived_const_cast
#undef derived_cast

#endif
