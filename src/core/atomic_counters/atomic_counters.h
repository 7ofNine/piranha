/***************************************************************************
 *   Copyright (C) 2009 by Francesco Biscani   *
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

#ifndef PIRANHA_ATOMIC_COUNTERS_H
#define PIRANHA_ATOMIC_COUNTERS_H

#include <cstddef>

#include "atomic_counter_msvc_long.h"

namespace piranha
{
        // TODO: make atomic counter for msvc generic using casting to long, see http://www.dlugosz.com/Repertoire/refman/Classics/atomic_counter_whitepaper.html
        // MSVC's atomic builtins operate on 32bit and 64bit with different naming conventions. But
        // if _WIN64 is defined the standard (non *64) methods also handle 64 bit parameters. The *64 methods are used
        // if _WIN64 is not defined but one wants to operate on 64 bit values. The MS documentation is competely misleading. See WinBase.h

        typedef atomic_counter_msvc_long atomic_counter_size_t;
}

#endif
