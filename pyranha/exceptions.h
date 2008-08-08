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

#ifndef PYRANHA_EXCEPTIONS_H
#define PYRANHA_EXCEPTIONS_H

#include <boost/python/exception_translator.hpp>

#include "../src/core/exceptions.h"

namespace pyranha
{
	template <class Exception>
	inline void exception_translator(const Exception &e)
	{
		PyErr_SetString(PyExc_UserWarning, e.what().c_str());
	}

	template <class Exception>
	inline void register_exception()
	{
		boost::python::register_exception_translator<Exception>(exception_translator<Exception>);
	}

	inline void translate_exceptions()
	{
		register_exception<piranha::bad_input>();
		register_exception<piranha::not_existing>();
		register_exception<piranha::not_implemented>();
		register_exception<piranha::unsuitable>();
		register_exception<piranha::division_by_zero>();
		register_exception<piranha::out_of_memory>();
	}
}

#endif
