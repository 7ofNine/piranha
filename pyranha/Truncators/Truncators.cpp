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


#include "../../src/core/truncators/degree.h"
#include "../../src/core/truncators/norm.h"
#include "../../src/core/truncators/truncators.h"
#include "../../src/core/mp.h"
#include "../commons.h"
#include "../exceptions.h"

#include "pybind11/pybind11.h"

#include <vector>
#include <string>


// Instantiate the pyranha Truncators module.
PYBIND11_MODULE(_Truncators, mt)
{
    //docstring_options docOptions(true, false, false);  //TODO: tob done
    pyranha::translate_exceptions();

    typedef void (*deg_set)(const int);
    typedef void (*q_deg_set)(const piranha::mp_rational &);
    typedef void (*p_deg_set)(const std::vector<std::string> &, const int);
    typedef void (*p_q_deg_set)(const std::vector<std::string> &, const piranha::mp_rational &);
    typedef void (*s_p_deg_set)(const std::string &, const int);
    typedef void (*s_p_q_deg_set)(const std::string &, const piranha::mp_rational &);


    pybind11::class_<piranha::truncators::Degree> td(mt, "__degree_truncator", "Minimum degree truncator.");
    td.def(pybind11::init<>());
    td.def("__repr__", &pyranha::py_print_to_string<piranha::truncators::Degree>);
    td.def_static("set", deg_set(&piranha::truncators::Degree::set), "Set truncation level of series' minimum degree to arg1.");
    td.def_static("set", s_p_deg_set(&piranha::truncators::Degree::set), "Set truncation level of series' partial minimum degree to arg2, "
        "relatively to Psym named arg1.");
    td.def_static("set", p_deg_set(&piranha::truncators::Degree::set), "Set truncation level of series' partial minimum degree to arg2, "
        "relatively to list of Psym names arg1.");
    td.def_static("set", q_deg_set(&piranha::truncators::Degree::set), "Set truncation level of series' minimum degree to arg1.");
    td.def_static("set", s_p_q_deg_set(&piranha::truncators::Degree::set), "Set truncation level of series' partial minimum degree to arg2, "
        "relatively to Psym named arg1.");
    td.def_static("set", p_q_deg_set(&piranha::truncators::Degree::set), "Set truncation level of series' partial minimum degree to arg2, "
        "relatively to list of Psym names arg1.");
    td.def_static("unset", &piranha::truncators::Degree::unset, "Clear minimum degree limit.");


    pybind11::class_<piranha::truncators::Norm> tn(mt, "__norm_truncator", "Norm truncator.");
    tn.def(pybind11::init<>());
    tn.def("__repr__", &pyranha::py_print_to_string<piranha::truncators::Norm>);
    tn.def_static("set", &piranha::truncators::Norm::set, "Set norm truncation level to arg1.");
    tn.def_static("unset", &piranha::truncators::Norm::unset, "Disable norm-based truncation.");

    mt.def("unset", &piranha::truncators::unset, "Unset all truncators.");
}
