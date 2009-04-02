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

#ifndef PIRANHA_INTEGER_TYPEDEFS_H
#define PIRANHA_INTEGER_TYPEDEFS_H

#include <boost/cstdint.hpp>

#include "config.h"

namespace piranha
{
	/// Integer selector.
	/**
	 * It detects wehther the platform is 64bit or 32bit, and sets maximum and minimum
	 * "fast" integer types accordingly, using the boost integer libraries. If the platform is
	 * other than 32bit or 64bit it won't define any type.
	 */
	template <int SizeOfPointer>
	struct int_selector {
		p_static_check(SizeOfPointer == 4, "Unknown pointer size.");
		typedef boost::int32_t  max_fast_int;
	};

	// Specialization for 64bit archs.
	template <>
	struct int_selector<8> {
		typedef boost::int64_t  max_fast_int;
	};

	/// Maximum fast integer, detected through piranha::int_selector.
	typedef int_selector < sizeof(void *) >::max_fast_int max_fast_int;
}

#endif
