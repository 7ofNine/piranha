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

#ifndef PIRANHA_CODED_MULTIPLIER_H
#define PIRANHA_CODED_MULTIPLIER_H

#include <boost/integer_traits.hpp>
#include <cstddef>

#include "../integer_typedefs.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	/// Toolbox for coded series multiplication.
	/**
	 * Intended to be inherited together with piranha::base_series_multiplier. It adds common methods for
	 * series multiplication through Kronecker codification.
	 */
	template <class Derived, class Key>
	class coded_series_multiplier
	{
			typedef max_fast_int coded_key_type;
			typedef boost::integer_traits<coded_key_type> traits;
			typedef typename Key::value_type value_type;
		public:
		protected:
			/// Is coded representation viable?
			bool							m_cr_is_viable;
			// Size of the coding vector, min_max vectors, etc.
			const std::size_t					m_size;
			// Vectors of minimum and maximum value pairs for the series being multiplied.
			std::vector<std::pair<int, int> >			m_min_max1;
			std::vector<std::pair<int, int> >			m_min_max2;
			// Vector of minimum and maximum value pairs for the resulting series.
			// GMP is used to avoid trespassing the range limits of max_fast_int.
			std::vector<std::pair<mp_integer, mp_integer> >		m_res_min_max;
			// Version of the above downcast to int.
			std::vector<std::pair<int, int> >			m_fast_res_min_max;
			// Coding vector.
			std::vector<max_fast_int>				m_coding_vector;
			// Mininum and maximum values of codes, and total number of codes.
			max_fast_int						m_h_min;
			max_fast_int						m_h_max;
			max_fast_int						m_h_tot;
			// Coded keys.
			std::vector<max_fast_int>				m_ckeys1;
			std::vector<max_fast_int>				m_ckeys2;
			// Densities of input series wrt the resulting series' codification.
			double							m_density1;
			double							m_density2;
	};
}

#undef derived_const_cast
#undef derived_cast

#endif
