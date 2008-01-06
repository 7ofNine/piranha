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

#include "config.h"

#define hard_assert(result) \
if (unlikely((result)==false)) \
{ \
  std::cout << __FILE__ << ':' << __LINE__ << " Assert failed" << std::endl; \
  exit(1); \
}

#if defined _ENABLE_ASSERTS

#define p_assert(result) hard_assert(result)

#define action_assert(action) \
if (unlikely((action)==false)) \
{ \
  std::cout << __FILE__ << ':' << __LINE__ << " Assert failed" << std::endl; \
  exit(1); \
}

#else

#define p_assert(__arg,...)

#define action_assert(action,...) action
#endif                                            // _ENABLE_ASSERTS
#endif
