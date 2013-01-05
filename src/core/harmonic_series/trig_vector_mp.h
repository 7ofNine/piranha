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

#ifndef PIRANHA_TRIG_VECTOR_MP_H
#define PIRANHA_TRIG_VECTOR_MP_H

#include <iostream>

#include "../mp.h"

namespace piranha
{
	inline void trig_vector_print_element_tex(std::ostream &out_stream, const mp_rational &q)
	{
		q.print_tex(out_stream);
	}


	template <class T>
	inline void trig_vector_print_element_tex(std::ostream &out_stream, const T &x)
	{
		out_stream << x;
	}


	inline double trig_vector_eval_element(const mp_rational &q)
	{
		return q.to_double();
	}


	template <class T>
	inline double trig_vector_eval_element(const T &x)
	{
		return static_cast<double>(x);
	}
}

#endif
