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

#ifndef PIRANHA_P_ASSERT_H
#define PIRANHA_P_ASSERT_H

#include <iostream>
#include <sstream>

#include "config.h" // For "unlikely()".
#include "exceptions.h"

#if defined _PIRANHA_ENABLE_ASSERTS
#define p_assert(result) \
	if (unlikely(!(result))) \
	{ \
		std::ostringstream tmp;\
		tmp << "Assertion failed at: " << __FILE__ << ", " << __LINE__ << ".\n"; \
		throw piranha::assertion_failure(tmp.str()); \
	}
#else
#define p_assert(__arg)
#endif // _PIRANHA_ENABLE_ASSERTS

#endif
