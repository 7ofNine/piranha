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

#ifndef PIRANHA_CONFIG_H
#define PIRANHA_CONFIG_H

#include <cmath>
#include <cstdlib>

#ifdef __GNUC__

#define GCC_VERSION (__GNUC__ * 100000 \
 + __GNUC_MINOR__ * 1000 \
 + __GNUC_PATCHLEVEL__ * 10)

#if GCC_VERSION < 304000
#error "The minimum GCC version required is 3.4"
#endif

#define likely(exp)   __builtin_expect(exp,1)
#define unlikely(exp) __builtin_expect(exp,0)

#else // __GNUC__

#warning Only the GNU compiler is officially supported.

#define likely(exp)   exp
#define unlikely(exp) exp

#endif // __GNUC__

// Platform switches.
#ifdef _PIRANHA_WIN32
  #define __ISNAN(x) _isnan(x)
  #define __JNL(n,x) jn(n,x)
  #define __ALIGNED_MALLOC(p,a,s) p=malloc(s)
  #define __PIRANHA_VISIBLE __declspec(dllexport)
#else
  #define __ISNAN(x) isnan(x)
  #define __JNL(n,x) jnl(n,x)
  #define __ALIGNED_MALLOC(p,a,s) posix_memalign(p,a,s)
  #define __PIRANHA_VISIBLE __attribute__ ((visibility("default")))
#endif

#endif
