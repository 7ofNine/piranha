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
#include <vector>
#include <string>

#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/module.hpp>
#include <boost/python/docstring_options.hpp>


#include "../../src/core/truncators/degree.h"
#include "../../src/core/truncators/norm.h"
#include "../../src/core/truncators/truncators.h"
#include "../../src/core/mp.h"
#include "../commons.h"
#include "../exceptions.h"

using namespace boost::python;
using namespace piranha;
using namespace pyranha;

// Instantiate the pyranha Truncators module.
BOOST_PYTHON_MODULE(_Truncators)
{
    docstring_options docOptions(true, false, false);
    translate_exceptions();

    typedef void (*deg_set)(const int);
    typedef void (*q_deg_set)(const mp_rational &);
    typedef void (*p_deg_set)(const std::vector<std::string> &, const int);
    typedef void (*p_q_deg_set)(const std::vector<std::string> &, const mp_rational &);
    typedef void (*s_p_deg_set)(const std::string &, const int);
    typedef void (*s_p_q_deg_set)(const std::string &, const mp_rational &);
    class_<truncators::Degree>("__degree_truncator", "Minimum degree truncator.", init<>())
    .def("__repr__", &py_print_to_string<truncators::Degree>)
    .def("set", deg_set(&truncators::Degree::set), "Set truncation level of series' minimum degree to arg1.")
    .def("set", s_p_deg_set(&truncators::Degree::set), "Set truncation level of series' partial minimum degree to arg2, "
        "relatively to Psym named arg1.")
    .def("set", p_deg_set(&truncators::Degree::set), "Set truncation level of series' partial minimum degree to arg2, "
        "relatively to list of Psym names arg1.")
    .def("set", q_deg_set(&truncators::Degree::set), "Set truncation level of series' minimum degree to arg1.")
    .def("set", s_p_q_deg_set(&truncators::Degree::set), "Set truncation level of series' partial minimum degree to arg2, "
        "relatively to Psym named arg1.")
    .def("set", p_q_deg_set(&truncators::Degree::set), "Set truncation level of series' partial minimum degree to arg2, "
        "relatively to list of Psym names arg1.").staticmethod("set")
    .def("unset", &truncators::Degree::unset, "Clear minimum degree limit.").staticmethod("unset");

    class_<truncators::Norm>("__norm_truncator", "Norm truncator.", init<>())
    .def("__repr__", &py_print_to_string<truncators::Norm>)
    .def("set", &truncators::Norm::set, "Set norm truncation level to arg1.").staticmethod("set")
    .def("unset", &truncators::Norm::unset, "Disable norm-based truncation.").staticmethod("unset");

    def("unset",&truncators::unset,"Unset all truncators.");
}
