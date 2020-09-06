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


#include "../../src/core/base_classes/named_series_def.h"
#include "../../src/core/config.h"
#include "../../src/core/mp.h"
#include "../../src/core/Psym.h"
#include "../../src/core/settings.h"
#include "../../src/core/stats.h"
#include "../args_tuple.h"
#include "../commons.h"
#include "../exceptions.h"
#include "../mp_classes.h"


#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include <cstddef>
#include <sstream>
#include <string>
#include <vector>
#include <functional>

std::string static inline py_psym_repr(const piranha::Psym &p)    //TODO: 
{                                                        // we should be able to use unicode!!!! e.g. for greek alphabet etc.       
    std::ostringstream stream;
    stream << "Symbol: '" << p.getName() << "' - " << "o:" << p.order() << " - e:[";
    const std::size_t size = p.getTimeEval().size();

    for (std::size_t i = 0; i < size; ++i) {
        stream << p.getTimeEval()[i];
        if (i < size - 1)
        {
            stream << ',';
        }
    }
    stream << ']';
    return stream.str();
}


static inline std::size_t py_psym_hash(const  piranha::Psym &p)
{
    return std::hash<std::string>()(p.getName());
}


static inline void  ed_set_item(piranha::EvalDict &d, const std::string &n, const double &value)
{
    d[n] = value;
}


static inline std::string py_stats_dump(const  piranha::stats &)
{
    return  piranha::stats::dump();
}



// Instantiate the pyranha Core module.
//BOOST_PYTHON_MODULE(_Core)
PYBIND11_MODULE(_Core, mc)
{   
//docstring_options docOptions(true, false, false);
pyranha::translate_exceptions();

// Interop between vectors of some types and Python tuples/lists.
//to_tuple_mapping<std::vector<std::string> >();

//from_python_sequence<std::vector<std::string>,variable_capacity_policy>();

//to_tuple_mapping<VectorPsym>();

//from_python_sequence<VectorPsym,variable_capacity_policy>();

//to_tuple_mapping<std::vector<VectorPsym> >();

//to_tuple_mapping<std::vector<double> >();

//from_python_sequence<std::vector<double>,variable_capacity_policy>();

//to_tuple_mapping<std::vector<mp_rational> >();

//from_python_sequence<std::vector<mp_rational>,variable_capacity_policy>();

//to_tuple_mapping<std::vector<mp_integer> >();

//from_python_sequence<std::vector<mp_integer>,variable_capacity_policy>();


// Expose evaluation dictionary.     TODO: disbaled until we need it
    pybind11::class_< piranha::EvalDict> ed(mc, "evaldict", "Evaluation dictionary");

    ed.def(pybind11::init<>());
    ed.def("__setitem__", &ed_set_item);    //TODO: is that all we need??? 

// Expose arguments tuples.   //is that used anywhere ????
    pyranha::expose_argsTuples<__PIRANHA_MAX_ECHELON_LEVEL>(mc);


    // MP classes.
    // class mp_rational
    pybind11::class_< piranha::mp_rational> mpr(pyranha::expose_real_mp_class< piranha::mp_rational>(mc, "rational", "Multi-precision rational number."));
    mpr.def(pybind11::init<const int &, const int &>());
    mpr.def(pybind11::init<const  piranha::mp_integer &, const  piranha::mp_integer &>());
    mpr.def_property_readonly("num", &piranha::mp_rational::get_num, "Numerator of the rational number");
    mpr.def_property_readonly("den", &piranha::mp_rational::get_den, "Denominator of the ratonal Number");
    mpr.def("choose", &piranha::mp_rational::choose, "Binomial coefficient (choose function).");


    // class mp_integer
    pybind11::class_< piranha::mp_integer> mpz(pyranha::expose_real_mp_class< piranha::mp_integer>(mc, "integer", "Multi-precision integer number."));
    mpz.def(pybind11::init<const  piranha::mp_rational &>());
   
    mpz.def("factorial", &piranha::mp_integer::factorial, "Factorial.");
    mpz.def("choose",    &piranha::mp_integer::choose,    "Binomial coefficient (choose function).");
    mpz.def("lcm",       &piranha::mp_integer::lcm,       "Set self to the least common multiplier of input arguments.");   //TODO: maybe a better interface. not just k.lcm(i,i) and for plain int i,j?
    mpz.def(pybind11::self %= piranha::mp_integer());                                                                       // todo: maybe also a better interface. make it assignable? only i%=j --> i works    
    mpz.def(pybind11::self %= int());


    // class settings (for pyranha)
    auto getMemoryLimit = [](pybind11::object) {return  piranha::settings::get_memory_limit(); };
    auto setMemoryLimit = [](pybind11::object, std::size_t v) {  piranha::settings::set_memory_limit(v); };
    auto getMemoryUsed  = [](pybind11::object) { return  piranha::settings::get_used_memory(); };
    auto getNthread     = [](pybind11::object) {return  piranha::settings::getNthread(); };
    auto setNthread     = [](pybind11::object, int n) { piranha::settings::setNthread(n); };
    auto getDebug       = [](pybind11::object) { return  piranha::settings::getDebug(); };
    auto setDebug       = [](pybind11::object, bool d) {  piranha::settings::setDebug(d); };
    auto getAlgorithm   = [](pybind11::object) { return  piranha::settings::getMultiplicationAlgorithm(); };
    auto setAlgorithm   = [](pybind11::object, piranha::settings::MultiplicationAlgorithm a) {  piranha::settings::setMultiplicationAlgorithm(a); };
   
    pybind11::class_< piranha::settings> set(mc, "__settings", "pyranha settings");
    set.def(pybind11::init<>());
    set.def_property_static("debug",                        getDebug, setDebug);  
    set.def_property_readonly_static("memory_used",         getMemoryUsed);                                                                         // , "amount of used memory in bytes.");                        // todo: how to do property with the docstring?????
    set.def_property_static("memory_limit",                 getMemoryLimit, setMemoryLimit);                                                        // todo: how to do property with the docstring?????
    set.def_property_static("max_pretty_print_size",        &piranha::settings::get_max_pretty_print_size, &piranha::settings::set_max_pretty_print_size);            // todo: how to do property with the docstring?????
    set.def_property_static("nthread",                      getNthread, setNthread);                                                                // todo: how to do property with the docstring?????
    set.def_property_static("multiplication_algorithm",     getAlgorithm, setAlgorithm);                                                            // todo: how to do property with the docstring?????
 

    // enum for multiplication algorithm   
    pybind11::enum_< piranha::settings::MultiplicationAlgorithm> enu(mc, "multiplication_algorithm", "Selected algorithm for multiplication");
    enu.value("automatic", piranha::settings::ALGORITHM_AUTOMATIC);
    enu.value("plain", piranha::settings::ALGORITHM_PLAIN);
    enu.value("vector_coded", piranha::settings::ALGORITHM_VECTOR_CODED);
    enu.value("hash_coded", piranha::settings::ALGORITHM_HASH_CODED);
    enu.export_values();


    // class psym
    // Should psym offer a latex method ?? not much really to show that way, probably just the name!
    pybind11::class_< piranha::Psym> ps(mc, "psym", "Symbol class.");
    
    ps.def(pybind11::init<const std::string&, const std::vector<double>&>());
    ps.def(pybind11::init<const std::string&, const double&>());
    ps.def(pybind11::init<const std::string&>());

    ps.def("__copy__", &pyranha::py_copy< piranha::Psym>);
    ps.def("__hash__", &py_psym_hash);
    ps.def("__repr__", &py_psym_repr);
    ps.def("eval", &piranha::Psym::eval);
    ps.def_property_readonly("name", &piranha::Psym::getName, pybind11::return_value_policy::copy);
    ps.def_property("time_eval", &piranha::Psym::getTimeEval, &piranha::Psym::setTimeEval, pybind11::return_value_policy::copy);
    ps.def_static("list", &piranha::Psym::list, "Get list of global psyms");
    ps.def_property("order", &piranha::Psym::order, &piranha::Psym::setOrder, "Get/set the order of the symbol for truncation (default: 1)");
    //
    ps.def(pybind11::self == pybind11::self);
    ps.def(pybind11::self != pybind11::self);


    // Stats class.
    pybind11::class_< piranha::stats> s(mc, "__stats", "Stats class.");
    
    s.def(pybind11::init<>());
    s.def("__repr__", &py_stats_dump);


}
