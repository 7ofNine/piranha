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

#ifndef PIRANHA_EXCEPTIONS_H
#define PIRANHA_EXCEPTIONS_H

#include <iostream>
#include <sstream>
#include <string>

#include "p_exceptions.h"

#include "config.h" // For "unlikely()".

#define PIRANHA_THROW(ex,s) P_EX_THROW(ex,s)

#if defined _PIRANHA_ENABLE_ASSERTS
#define PIRANHA_ASSERT(result) \
	if (unlikely(!(result))) \
	{ \
		std::ostringstream tmp;\
		tmp << "Assertion failed at: " << __FILE__ << ", " << __LINE__ << ".\n"; \
		tmp << "Assertion failures should not happen and denote a bug.\n"; \
		tmp << "Please help us correct the problem by sending a bug report to piranha-psm@googlegroups.com\n"; \
		tmp << "providing as many details as possible. Thanks!\n"; \
		PIRANHA_THROW(assertion_error,"assertion error"); \
	}
#else
#define PIRANHA_ASSERT(__arg)
#endif // _PIRANHA_ENABLE_ASSERTS

#endif
