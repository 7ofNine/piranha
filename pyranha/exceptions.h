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

#include "../src/core/p_exceptions.h"


#include "pybind11/pybind11.h"

namespace pyranha
{   

    inline void translate_exceptions()
    {

        pybind11::register_exception_translator([](std::exception_ptr p)
        {
            
            try {
                if (p) std::rethrow_exception(p);
            }
            catch (const index_error &e){
                PyErr_SetString(PyExc_IndexError, e.what());
            }
            catch(const value_error &e){

                PyErr_SetString(PyExc_ValueError, e.what());
            }
            catch(const type_error &e)
            {
                PyErr_SetString(PyExc_TypeError, e.what());
            }
            catch(const assertion_error &e)
            {
                PyErr_SetString(PyExc_AssertionError, e.what());
            }
            catch(const not_implemented_error &e)
            {
                PyErr_SetString(PyExc_NotImplementedError, e.what());
            }
            catch(const memory_error &e)
            {
                PyErr_SetString(PyExc_MemoryError, e.what());
            }
            catch(const zero_division_error &e)
            {
                PyErr_SetString(PyExc_ZeroDivisionError, e.what());
            }
            catch(const boost::numeric::bad_numeric_cast &e)
            {
                PyErr_SetString(PyExc_OverflowError, e.what());
            }

            });

    }
}

#endif
