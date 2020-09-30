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

#ifndef PIRANHA_NUMERICAL_CONTAINER_MP_H
#define PIRANHA_NUMERICAL_CONTAINER_MP_H

#include <boost/lexical_cast.hpp>

#include "numerical_container_tag.h"

#include <complex>
#include <iostream>
#include <string>
#include <type_traits>



namespace piranha
{
	template <class T>
	struct numerical_container_eval_type_determiner
	{
		typedef double type;
	};

	template <class T>
	struct numerical_container_eval_type_determiner<std::complex<T> >
	{
		typedef std::complex<double> type;
	};

	template <typename T>
	struct numerical_container_constructor_selector
	{
		static const T &run(const T &x)
		{
			return x;
		}
	};

	template <class T> requires std::is_base_of_v<numerical_container_tag, T>
	struct numerical_container_constructor_selector<T> 
	{
		static const typename T::numerical_type &run(const T &x)
		{
			return x.get_value();
		}
	};

	template <class Value>
	struct numerical_container_print_pretty_selector
	{
		static void run(const Value &x, std::ostream &os)
		{
			os << boost::lexical_cast<std::string>(x);
		}
	};

	template <class RealValue>
	struct numerical_container_print_pretty_selector<std::complex<RealValue> >
	{
		static void run(const std::complex<RealValue> &c, std::ostream &os)
		{
			os << '(';
			numerical_container_print_pretty_selector<RealValue>::run(c.real(),os);
			if (c.imag() >= 0) {
				os << '+';
			}
			numerical_container_print_pretty_selector<RealValue>::run(c.imag(),os);
			os << "j)";
		}
	};

	template <class Value>
	struct numerical_container_print_tex_selector
	{
		static void run(const Value &x, std::ostream &os)
		{
			os << ' ';
			os << boost::lexical_cast<std::string>(x);
			os << ' ';
		}
	};

	template <class RealValue>
	struct numerical_container_print_tex_selector<std::complex<RealValue> >
	{
		static void run(const std::complex<RealValue> &c, std::ostream &os)
		{
			os << "\\left(";
			numerical_container_print_tex_selector<RealValue>::run(c.real(),os);
			if (c.imag() >= 0) {
				os << '+';
			}
			numerical_container_print_tex_selector<RealValue>::run(c.imag(),os);
			os << "\\, \\imath\\right)";
		}
	};

	// Automatically handle complex - scalar comparisons.
	template <class T, class U>
	struct numerical_container_value_comparison
	{
		static bool run(const T &x, const U &y)
		{
			return (x == y);
		}
	};

	template <class T, class U>
	struct numerical_container_value_comparison<std::complex<T>,U>
	{
		static bool run(const std::complex<T> &c, const U &x)
		{
			return (c.real() == x && c.imag() == 0);
		}
	};


	template <typename Derived, typename T>
	struct numerical_container_equality_selector
	{
		static bool run(const Derived &cf, const T &x)
		{
			if constexpr (std::is_base_of_v<numerical_container_tag, T>)
			{
				return numerical_container_value_comparison<typename Derived::numerical_type, typename T::numerical_type>::run(cf.m_value, x.get_value());
			}
			else
			{
				return numerical_container_value_comparison<typename Derived::numerical_type, T>::run(cf.m_value, x);
			}
		}
	};


	template <typename T>
	struct numerical_container_add_selector
	{
		template <typename Derived>
		static Derived &run(Derived &cf, const T &x)
		{
			if constexpr (std::is_base_of_v<numerical_container_tag, T>)
			{
				cf.m_value += x.get_value();
			}
			else
			{
				cf.m_value += x;
			}
			
			return cf;
		}
	};


	template <typename T>
	struct numerical_container_subtract_selector
	{
		template <typename Derived>
		static Derived& run(Derived& cf, const T& x)
		{
			if constexpr (std::is_base_of_v<numerical_container_tag, T>)
			{
				cf.m_value -= x.get_value();
			}
			else
			{
				cf.m_value -= x;
			}
			return cf;
		}
	};


	template <typename T, typename U>
	struct numerical_container_value_multiplication
	{
		static void run(T &x, const U &y)
		{
			x *= y;
		}
	};


	template <typename T>
	struct numerical_container_multiply_selector
	{
		template <typename Derived>
		static Derived &run(Derived &cf, const T &x)
		{
			if constexpr (std::is_base_of_v<numerical_container_tag, T>)
			{
				numerical_container_value_multiplication<typename Derived::numerical_type,typename T::numerical_type>::run(cf.m_value, x.get_value());
			}
			else
			{
				numerical_container_value_multiplication<typename Derived::numerical_type, T>::run(cf.m_value, x);
			}
			return cf;
		}
	};


	template <typename T, typename U>
	struct numerical_container_value_division
	{
		static void run(T &x, const U &y)
		{
			x /= y;
		}
	};


	template <typename T>
	struct numerical_container_divide_selector
	{
		template <typename Derived>
		static Derived& run(Derived& cf, const T& x)
		{
			if constexpr ( std::is_base_of_v<numerical_container_tag, T>)
			{
				numerical_container_value_division<typename Derived::numerical_type, typename T::numerical_type>::run(cf.m_value, x.get_value());
			}
			else
			{
				numerical_container_value_division<typename Derived::numerical_type, T>::run(cf.m_value, x);
			}
			return cf;
		}
	};

}

#endif
