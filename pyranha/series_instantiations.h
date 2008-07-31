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

#ifndef PYRANHA_SERIES_INSTANTIATIONS_H
#define PYRANHA_SERIES_INSTANTIATIONS_H

#include <boost/lexical_cast.hpp>
#include <boost/python/class.hpp>
#include <boost/python/iterator.hpp>
#include <boost/python/operators.hpp>
#include <boost/python/pure_virtual.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <complex>
#include <string>
#include <utility> // For std::pair.

#include "../src/core/integer_typedefs.h"
#include "../src/core/psym.h"
#include "../src/core/type_traits.h"
#include "cf_key_bindings.h"

namespace pyranha
{
	template <class Series>
	inline typename Series::const_iterator py_series_begin(const Series &s)
	{
		return s.begin();
	}

	template <class Series>
	inline typename Series::const_iterator py_series_end(const Series &s) {
		return s.end();
	}

	template <class Series>
	inline Series py_series_copy(const Series &s)
	{
		return Series(s);
	}

	template <class Series, class Term>
	inline void py_series_append(Series &s, const Term &t) {
		s.insert(t,s.arguments());
	}

	/// Basic series instantiation.
	template <class T>
	std::pair<boost::python::class_<T>, boost::python::class_<typename T::term_type> >
	series_basic_instantiation(const std::string &name, const std::string &description)
	{
		typedef typename T::term_type term_type;
		// Expose the term type.
		boost::python::class_<term_type> term_inst((name + "_term").c_str(),
				(std::string("Term for: ") + description).c_str());
		term_inst.def_readonly("cf", &term_type::m_cf);
		term_inst.def_readonly("key", &term_type::m_key);
		// Expose the manipulator class.
		boost::python::class_<T> inst(name.c_str(), description.c_str());
		inst.def(boost::python::init<const T &>());
		inst.def(boost::python::init<const std::string &>());
		inst.def(boost::python::init<const piranha::max_fast_int &>());
		inst.def(boost::python::init<const double &>());
		inst.def(boost::python::init<const piranha::psym &>());
		// Some special methods.
		inst.def("__copy__", &py_series_copy<T>);
		inst.def("__iter__", boost::python::range(&py_series_begin<T>, &py_series_end<T>));
		inst.def("__len__", &T::length);
		inst.def("__repr__", &T::py_repr);
		// Pyranha-specific special methods.
		inst.add_property("__arguments_description__", &T::py_arguments_description);
		inst.add_property("__arguments__", &T::py_arguments);
		inst.def("__set_arguments__", &T::set_arguments);
		inst.def("__set_shared_arguments__", &T::py_shared_arguments_set);
		typedef void (T::*trim_free)();
		inst.def("__trim__", trim_free(&T::trim));
		inst.def("__append__", &py_series_append<T,term_type>);
		inst.def("save_to", &T::save_to, "Save series to file.");
		typedef typename T::eval_type (T::*eval_free)(const double &) const;
		inst.def("eval", eval_free(&T::eval));
		typedef double (T::*norm_named)() const;
		inst.def("norm", norm_named(&T::norm));
		inst.def("atoms", &T::atoms);
		inst.def("swap", &T::swap);
		// NOTICE: the order seems important here, if we place *=int before *=double we
		// will get just *=double in Python. Go figure...
		// Addition and subtraction.
		inst.def(boost::python::self += piranha::max_fast_int());
		inst.def(boost::python::self += double());
		inst.def(boost::python::self += boost::python::self);
		inst.def(boost::python::self + piranha::max_fast_int());
		inst.def(piranha::max_fast_int() + boost::python::self);
		inst.def(boost::python::self + double());
		inst.def(double() + boost::python::self);
		inst.def(boost::python::self + boost::python::self);
		inst.def(boost::python::self -= piranha::max_fast_int());
		inst.def(boost::python::self -= double());
		inst.def(boost::python::self -= boost::python::self);
		inst.def(boost::python::self - piranha::max_fast_int());
		inst.def(piranha::max_fast_int() - boost::python::self);
		inst.def(boost::python::self - double());
		inst.def(double() - boost::python::self);
		inst.def(boost::python::self - boost::python::self);
		inst.def(-boost::python::self);
		// Multiplication.
		inst.def(boost::python::self *= piranha::max_fast_int());
		inst.def(boost::python::self *= double());
		inst.def(boost::python::self *= boost::python::self);
		inst.def(boost::python::self*piranha::max_fast_int());
		inst.def(piranha::max_fast_int()*boost::python::self);
		inst.def(boost::python::self*double());
		inst.def(double()*boost::python::self);
		inst.def(boost::python::self*boost::python::self);
		// Division.
		inst.def(boost::python::self /= piranha::max_fast_int());
		inst.def(boost::python::self /= double());
		inst.def(boost::python::self / piranha::max_fast_int());
		inst.def(boost::python::self / double());
		// Exponentiation.
		typedef T(T::*pow_double)(const double &) const;
		typedef T(T::*pow_int)(const piranha::max_fast_int &) const;
		inst.def("__pow__", pow_double(&T::pow));
		inst.def("__pow__", pow_int(&T::pow));
		typedef T(T::*named_root)(const piranha::max_fast_int &) const;
		inst.def("root", named_root(&T::root));
		return std::make_pair(inst, term_inst);
	}

	template <class T>
	void series_complex_instantiation(boost::python::class_<std::complex<T> > &instc, boost::python::class_<T> &inst)
	{
		// Ctors.
		instc.def(boost::python::init<const std::complex<piranha::max_fast_int> &>());
		instc.def(boost::python::init<const std::complex<double> &>());
		instc.def(boost::python::init<const T &>());
		instc.def(boost::python::init<const T &, const T &>());
		inst.def("complex", &T::complex, "Return complex counterpart.");
		// Addition and subtraction.
		instc.def(boost::python::self += std::complex<piranha::max_fast_int>());
		instc.def(boost::python::self += std::complex<double>());
		instc.def(boost::python::self += T());
		instc.def(boost::python::self + piranha::max_fast_int());
		instc.def(std::complex<piranha::max_fast_int>() + boost::python::self);
		instc.def(boost::python::self + std::complex<double>());
		instc.def(std::complex<double>() + boost::python::self);
		instc.def(boost::python::self + T());
		instc.def(T() + boost::python::self);
		instc.def(boost::python::self -= std::complex<piranha::max_fast_int>());
		instc.def(boost::python::self -= std::complex<double>());
		instc.def(boost::python::self -= T());
		instc.def(boost::python::self - piranha::max_fast_int());
		instc.def(std::complex<piranha::max_fast_int>() - boost::python::self);
		instc.def(boost::python::self - std::complex<double>());
		instc.def(std::complex<double>() - boost::python::self);
		instc.def(boost::python::self - T());
		instc.def(T() - boost::python::self);
		// Multiplication.
		instc.def(boost::python::self *= std::complex<piranha::max_fast_int>());
		instc.def(boost::python::self *= std::complex<double>());
		instc.def(boost::python::self *= T());
		instc.def(boost::python::self * piranha::max_fast_int());
		instc.def(std::complex<piranha::max_fast_int>() * boost::python::self);
		instc.def(boost::python::self * std::complex<double>());
		instc.def(std::complex<double>() * boost::python::self);
		instc.def(boost::python::self * T());
		instc.def(T() * boost::python::self);
		// Division.
		instc.def(boost::python::self /= std::complex<piranha::max_fast_int>());
		instc.def(boost::python::self /= std::complex<double>());
		instc.def(boost::python::self / std::complex<piranha::max_fast_int>());
		instc.def(boost::python::self / std::complex<double>());
		// Real and imaginary parts assignment and extraction.
		typedef T(std::complex<T>::*comp_get)() const;
		typedef std::complex<T> &(std::complex<T>::*comp_set)(const T &);
		instc.def("real", comp_get(&std::complex<T>::real), "Get real part.");
		instc.def("imag", comp_get(&std::complex<T>::imag), "Get imaginary part.");
		instc.def("real", comp_set(&std::complex<T>::real),
				  boost::python::return_internal_reference<1>(), "Set real part.");
		instc.def("imag", comp_set(&std::complex<T>::imag),
				  boost::python::return_internal_reference<1>(), "Set imaginary part.");
	}

	template <class T>
	void series_trigonometric_instantiation(boost::python::class_<T> &inst)
	{
		typedef std::complex<T> (T::*named_complexp)() const;
		inst.def("complexp", named_complexp(&T::complexp));
		inst.def("cos", &T::cos);
		inst.def("sin", &T::sin);
	}

	template <class T>
	T py_series_partial_name(const T &series, const std::string &p_name)
	{
		return series.partial(piranha::psyms::get(p_name));
	}

	template <class T>
	T py_series_partial_psym(const T &series, const piranha::psym &p)
	{
		return series.partial(p);
	}

	template <class T>
	void series_differential_instantiation(boost::python::class_<T> &inst)
	{
		inst.def("partial", &py_series_partial_psym<T>);
		inst.def("partial", &py_series_partial_name<T>);
	}

	template <class T, class Series>
	T py_series_sub_psym_psym(const T &series, const piranha::psym &p, const piranha::psym &q)
	{
		return series.template sub<Series>(p, Series(q));
	}

	template <class T, class Series>
	T py_series_sub_psym_series(const T &series, const piranha::psym &p, const Series &sub)
	{
		return series.template sub<Series>(p, sub);
	}

	template <class T, class Series>
	T py_series_sub_name_name(const T &series, const std::string &p_name, const std::string &q_name)
	{
		return series.template sub<Series>(piranha::psyms::get(p_name), Series(piranha::psyms::get(q_name)));
	}

	template <class T, class Series>
	T py_series_sub_name_series(const T &series, const std::string &p_name, const Series &sub)
	{
		return series.template sub<Series>(piranha::psyms::get(p_name), sub);
	}

	template <class T, class Series>
	void series_sub_instantiation(boost::python::class_<T> &inst)
	{
		inst.def("sub", py_series_sub_name_name<T, Series>);
		inst.def("sub", py_series_sub_psym_psym<T, Series>);
		inst.def("sub", py_series_sub_name_series<T, Series>);
		inst.def("sub", py_series_sub_psym_series<T, Series>);
	}

	template <class T>
	void series_special_functions_instantiation(boost::python::class_<T> &inst)
	{
		inst.def("besselJ", &T::besselJ, "Bessel function of the first kind of integer order.");
		inst.def("dbesselJ", &T::dbesselJ, "Partial derivative of Bessel function of the first kind of integer order.");
	}

#define __celmec_inst(arg) \
	inst.def(#arg,from_series(&T::arg),arg##_docstring) \
	.def(#arg,from_psym(&T::arg),arg##_docstring) \
	.def(#arg,from_name(&T::arg),arg##_docstring) \
	.staticmethod(#arg);

	template <class T>
	void celmec_instantiation(boost::python::class_<T> &inst)
	{
		typedef T(*from_series)(const T &, const T &);
		typedef T(*from_psym)(const piranha::psym &, const piranha::psym &);
		typedef T(*from_name)(const std::string &, const std::string &);
		const char *r_a_docstring = "Elliptic expansion of r / a.";
		const char *cos_f_docstring = "Elliptic expansion of cos(f).";
		const char *sin_f_docstring = "Elliptic expansion of sin(f).";
		const char *cos_E_docstring = "Elliptic expansion of cos(E).";
		const char *sin_E_docstring = "Elliptic expansion of sin(E).";
		const char *E_docstring = "Elliptic expansion of E (solution of Kepler's equation).";
		__celmec_inst(r_a);
		__celmec_inst(cos_f);
		__celmec_inst(sin_f);
		__celmec_inst(cos_E);
		__celmec_inst(sin_E);
		__celmec_inst(E);
	}

#undef __celmec_inst

	template <class T>
	void power_series_instantiation(boost::python::class_<T> &inst)
	{
		inst.def("degree", &T::degree, "Get the degree of the power series.");
		inst.def("min_degree", &T::min_degree, "Get the minimum degree of the power series.");
	}

	template <class T>
	void common_polynomial_instantiation(boost::python::class_<T> &inst)
	{
		series_differential_instantiation(inst);
		series_special_functions_instantiation(inst);
		power_series_instantiation(inst);
	}

	template <class T>
	void common_poisson_series_instantiation(boost::python::class_<T> &inst, const std::string &name)
	{
		series_differential_instantiation(inst);
		series_special_functions_instantiation(inst);
		power_series_instantiation(inst);
		// Expose the polynomial coefficient.
		typedef typename T::term_type::cf_type cf_type;
		cf_bindings<cf_type>((name + "_cf").c_str(), "")
		.def("degree", &cf_type::degree)
		.def("min_degree", &cf_type::min_degree);
	}

	template <class T>
	void common_fourier_series_instantiation(boost::python::class_<T> &inst)
	{
		series_special_functions_instantiation(inst);
		series_differential_instantiation(inst);
	}
}

#endif
