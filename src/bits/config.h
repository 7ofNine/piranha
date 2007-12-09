/***************************************************************************
 *   Copyright (C) 2007 by Francesco Biscani   *
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

#ifndef PIRANHA_CONFIG_H
#define PIRANHA_CONFIG_H

#ifndef __GNUC__
#error "GCC is the only supported compiler"
#endif

#define GCC_VERSION (__GNUC__ * 100000 \
  + __GNUC_MINOR__ * 1000 \
  + __GNUC_PATCHLEVEL__ * 10)

#if GCC_VERSION < 304000
#error "Minimum required GCC version is 3.4"
#endif

#if GCC_VERSION < 402000
#include "hash_set_hm.h"
#else
#include "unordered_set_hm.h"
#endif

#define likely(exp)   __builtin_expect(exp,1)
#define unlikely(exp) __builtin_expect(exp,0)

#define _PIRANHA_DISPLAY_PROGRESS (true)

#endif
