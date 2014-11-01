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

#ifndef PIRANHA_EXPO_VECTOR_MP_H
#define PIRANHA_EXPO_VECTOR_MP_H

#include <boost/numeric/conversion/cast.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/utility/enable_if.hpp>
#include <iostream>

#include "../exceptions.h"
#include "../mp.h"

namespace piranha
{
	inline void expoVectorPrintElementPretty(std::ostream &outStream, const mp_rational &q)
	{
		const bool needBracket = (q.get_den() != 1);
		if (needBracket)
        {
			outStream << '(';
		}
		outStream << q;

		if (needBracket)
        {
			outStream << ')';
		}
	}


	template <class T>
	inline void expoVectorPrintElementPretty(std::ostream &outStream, const T &x)
	{
		outStream << x;
	}

	inline void expoVectorPrintElementTEX(std::ostream &outStream, const mp_rational &q)
	{
		q.print_tex(outStream);
	}

	template <class T>
	inline void expoVectorPrintElementTEX(std::ostream &outStream, const T &x)
	{
		outStream << x;
	}

	// Default behaviour is to refuse to do anything is expo vector is not unity.
	template <class T, class Enable = void>
	struct ExpoVectorPowDouble
	{
		template <class ExpoVector>
		static ExpoVector run(const ExpoVector &expoVector, const double &)
		{
			if (!expoVector.isUnity())
            {
				PIRANHA_THROW(value_error, "cannot raise non-unity exponent vector to real power");
			}

			return expoVector;
		}
	};


	// For rationals, we build a rational from the double and go on with the exponentiation.
	template <class T>
	struct ExpoVectorPowDouble<T, typename boost::enable_if<boost::is_same<T, mp_rational> >::type>
	{
		template <class ExpoVector>
		static ExpoVector run(const ExpoVector &expoVector, const double &x)
		{
			return expoVector.genericPow(mp_rational(x));
		}
	};


	// Default behaviour is to refuse to do anything is resulting expo vector is not unity.
	template <class T, class Enable = void>
	struct ExpoVectorPowRational
	{
		template <class ExpoVector>
		static ExpoVector run(const ExpoVector &expoVector, const mp_rational &q)
		{
			typedef typename ExpoVector::size_type  size_type;
			typedef typename ExpoVector::value_type value_type;

			ExpoVector retval(expoVector);
			const size_type size = expoVector.size();

			for(size_type i = 0; i < size; ++i)
            {
				mp_rational tmp(q);
				tmp *= expoVector[i];
				if (tmp.get_den() != 1) 
                {
					PIRANHA_THROW(value_error, "exponent is not suitable for the calculation of rational power");
				}

				retval[i] = boost::numeric_cast<value_type>(tmp.to_int());
			}

			return retval;
		}
	};


	// For rationals, we do simple exponentiation.
	template <class T>
	struct ExpoVectorPowRational<T, typename boost::enable_if<boost::is_same<T, mp_rational> >::type>
	{
		template <class ExpoVector>
		static ExpoVector run(const ExpoVector &expoVector, const mp_rational &q)
		{
			return expoVector.genericPow(q);
		}
	};
}

#endif
