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
#include <utility>
#include <vector>

#include "../integer_typedefs.h"
#include "../mp.h"

#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)

namespace piranha
{
	/// Toolbox for coded series multiplication.
	/**
	 * Intended to be inherited together with piranha::base_series_multiplier. It adds common methods for
	 * series multiplication through Kronecker codification.
	 */
	template <class Derived, class Key, bool SubtractKeys>
	class coded_series_multiplier
	{
			typedef boost::integer_traits<max_fast_int> traits;
		public:
			/// Value type of the key.
			typedef typename Key::value_type value_type;
		protected:
			/// Is coded representation viable?
			bool							m_cr_is_viable;
			/// Size of the coding vector, min_max vectors, etc.
			const std::size_t					m_size;
			/// Vectors of minimum and maximum value pairs for the first series.
			std::vector<std::pair<value_type, value_type> >		m_min_max1;
			/// Vectors of minimum and maximum value pairs for the second series.
			std::vector<std::pair<value_type, value_type> >		m_min_max2;
			/// Vector of minimum and maximum value pairs for the resulting series, MP version.
			std::vector<std::pair<mp_integer, mp_integer> >		m_res_min_max;
			/// Vector of minimum and maximum value pairs for the resulting series, value_type version.
			std::vector<std::pair<value_type, value_type> >		m_fast_res_min_max;
			/// Coding vector.
			std::vector<max_fast_int>				m_coding_vector;
			/// Mininum code.
			max_fast_int						m_h_min;
			/// Maximum code.
			max_fast_int						m_h_max;
			/// Total number of codes.
			max_fast_int						m_h_tot;
			/// Coded keys for the first series.
			std::vector<max_fast_int>				m_ckeys1;
			/// Coded keys for the second series.
			std::vector<max_fast_int>				m_ckeys2;
			/// Density of the first series.
			double							m_density1;
			/// Density of the second series.
			double							m_density2;
	};
}

#undef derived_const_cast
#undef derived_cast

#endif
