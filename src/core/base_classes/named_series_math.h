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

#include "../exceptions.h"
#include "../integer_typedefs.h"
#include "../ntuple.h"
#include "../psym.h"
#include "../settings.h"

namespace piranha
{
	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	template <bool Sign, class Derived2>
	inline Derived &named_series<__PIRANHA_NAMED_SERIES_TP>::merge_with_series(const Derived2 &s2)
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
	inline Derived &named_series<__PIRANHA_NAMED_SERIES_TP>::mult_by_series(const Derived2 &s2)
	{
		// First we merge the arguments of the two series.
		merge_args(s2);
		// Then we perform the multiplication.
		derived_cast->mult_by(s2, m_arguments);
		trim();
		return *derived_cast;
	}

	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	template <bool Sign, class Number>
	inline Derived &named_series<__PIRANHA_NAMED_SERIES_TP>::merge_number_helper(const Number &x)
	{
		derived_cast->template merge_with_number<Sign>(x, m_arguments);
		trim();
		return *derived_cast;
	}

	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline Derived &named_series<__PIRANHA_NAMED_SERIES_TP>::operator+=(const max_fast_int &n)
	{
		return merge_number_helper<true>(n);
	}

	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline Derived &named_series<__PIRANHA_NAMED_SERIES_TP>::operator+=(const double &x)
	{
		return merge_number_helper<true>(x);
	}

	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline Derived &named_series<__PIRANHA_NAMED_SERIES_TP>::operator+=(const Derived &s2)
	{
		return merge_with_series<true>(s2);
	}

	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline Derived &named_series<__PIRANHA_NAMED_SERIES_TP>::operator-=(const max_fast_int &n)
	{
		return merge_number_helper<false>(n);
	}

	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline Derived &named_series<__PIRANHA_NAMED_SERIES_TP>::operator-=(const double &x)
	{
		return merge_number_helper<false>(x);
	}

	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline Derived &named_series<__PIRANHA_NAMED_SERIES_TP>::operator-=(const Derived &s2)
	{
		return merge_with_series<false>(s2);
	}

	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline Derived named_series<__PIRANHA_NAMED_SERIES_TP>::operator-() const
	{
		Derived retval(*derived_const_cast);
		retval *= (max_fast_int)(-1);
		return retval;
	}

	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	template <class Number>
	inline Derived &named_series<__PIRANHA_NAMED_SERIES_TP>::mult_number_helper(const Number &x)
	{
		derived_cast->mult_by(x, m_arguments);
		trim();
		return *derived_cast;
	}

	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline Derived &named_series<__PIRANHA_NAMED_SERIES_TP>::operator*=(const max_fast_int &n)
	{
		return mult_number_helper(n);
	}

	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline Derived &named_series<__PIRANHA_NAMED_SERIES_TP>::operator*=(const double &x)
	{
		return mult_number_helper(x);
	}

	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline Derived &named_series<__PIRANHA_NAMED_SERIES_TP>::operator*=(const Derived &s2)
	{
		return mult_by_series(s2);
	}

	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	template <class Number>
	inline Derived &named_series<__PIRANHA_NAMED_SERIES_TP>::divide_number_helper(const Number &x)
	{
		derived_cast->divide_by(x, m_arguments);
		trim();
		return *derived_cast;
	}

	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline Derived &named_series<__PIRANHA_NAMED_SERIES_TP>::operator/=(const max_fast_int &n)
	{
		return divide_number_helper(n);
	}

	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline Derived &named_series<__PIRANHA_NAMED_SERIES_TP>::operator/=(const double &x)
	{
		return divide_number_helper(x);
	}

	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline Derived named_series<__PIRANHA_NAMED_SERIES_TP>::choose(const max_fast_int &n, const max_fast_int &k)
	{
		return Derived::choose(n,k,args_tuple_type());
	}

	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline Derived named_series<__PIRANHA_NAMED_SERIES_TP>::factorial(const max_fast_int &n)
	{
		return Derived::factorial(n,args_tuple_type());
	}

	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline Derived named_series<__PIRANHA_NAMED_SERIES_TP>::pow(const double &x) const
	{
		Derived retval(derived_const_cast->pow(x, m_arguments));
		retval.m_arguments = m_arguments;
		retval.trim();
		return retval;
	}

	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline Derived named_series<__PIRANHA_NAMED_SERIES_TP>::pow(const max_fast_int &n) const
	{
		Derived retval(derived_const_cast->pow(n, m_arguments));
		retval.m_arguments = m_arguments;
		retval.trim();
		return retval;
	}

	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline Derived named_series<__PIRANHA_NAMED_SERIES_TP>::root(const max_fast_int &n) const
	{
		Derived retval(derived_const_cast->root(n, m_arguments));
		retval.m_arguments = m_arguments;
		retval.trim();
		return retval;
	}

	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline Derived named_series<__PIRANHA_NAMED_SERIES_TP>::inv() const
	{
		Derived retval = derived_const_cast->inv_(m_arguments);
		retval.m_arguments = m_arguments;
		retval.trim();
		return retval;
	}

	/// Partial derivative with respect to a piranha::psym.
	template <__PIRANHA_NAMED_SERIES_TP_DECL>
	inline Derived named_series<__PIRANHA_NAMED_SERIES_TP>::partial(const psym &arg, const max_fast_int &n) const
	{
		typedef typename ntuple<std::pair<bool, size_t>, n_arguments_sets>::type pos_tuple_type;
		pos_tuple_type pos_tuple;
		psym_p p = psyms::get_pointer(arg);
		named_series_get_psym_p_positions<pos_tuple_type, args_tuple_type>::run(p, pos_tuple, m_arguments);
		Derived retval(derived_const_cast->partial(n, pos_tuple, m_arguments));
		retval.m_arguments = m_arguments;
		retval.trim();
		return retval;
	}
}

#endif
