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

#ifndef PIRANHA_ATOMIC_COUNTER_H
#define PIRANHA_ATOMIC_COUNTER_H

#include "config.h"

#ifdef _PIRANHA_MT
#if defined( __GNUC__ ) && GCC_VERSION >= 401000

#include "atomic_counter_gcc_41.h"

namespace piranha
{
	typedef atomic_counter_gcc_41<size_t> unsigned_atomic_counter;
}

#else // Not GCC or GCC < 4.1.

#include "atomic_counter_generic.h"

namespace piranha
{
	typedef atomic_counter_generic<size_t> unsigned_atomic_counter;
}

#endif // Compiler selection in case of MT.
#else // _PIRANHA_MT

namespace piranha
{
	// If multi-thread support was not enabled use plain old integers as counters.
	typedef size_t unsigned_atomic_counter;
}

#endif // _PIRANHA_MT

#endif
