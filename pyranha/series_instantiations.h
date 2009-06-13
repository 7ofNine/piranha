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
#include <complex>
#include <string>
#include <utility> // For std::pair.

#include "../src/core/mp.h"
#include "../src/core/psym.h"
#include "../src/core/type_traits.h"
#include "cf_key_bindings.h"
#include "commons.h"

namespace pyranha
{
	/// Basic series instantiation.
	template <class T>
	inline std::pair<boost::python::class_<T>, boost::python::class_<typename T::term_type> >
	series_basic_instantiation(const std::string &name, const std::string &description)
	{
		typedef typename T::term_type term_type;
		// Expose the term type.
		boost::python::class_<term_type> term_inst(py_series_term<term_type>(name,description));
		// Expose the manipulator class.
		boost::python::class_<T> inst(name.c_str(), description.c_str());
		inst.def(boost::python::init<const T &>());
		inst.def(boost::python::init<const std::string &>());
		inst.def(boost::python::init<const double &>());
		inst.def(boost::python::init<const piranha::mp_rational &>());
		inst.def(boost::python::init<const piranha::mp_integer &>());
		inst.def(boost::python::init<const piranha::psym &>());
		// Some special methods.
		inst.def("__copy__", &py_copy<T>);
		inst.def("__iter__", boost::python::iterator<T>());
		inst.def("__len__", &T::length);
		inst.def("__repr__", &py_print_to_string<T>);
		// Pyranha-specific special methods.
		inst.def("__psi__", &psi0<T>);
		inst.def("__psi__", &psi1<T>);
		inst.def("__psi__", &psi2<T>);
		inst.add_property("__arguments_description__", &py_series_arguments_description<T>);
		inst.add_property("__arguments__", &py_series_arguments<T>);
		inst.def("__set_arguments__", &T::set_arguments);
		//typedef void (T::*trim_free)();
		//inst.def("__trim__", trim_free(&T::trim));
		//inst.def("__append__", &py_series_append<T,term_type>);
		inst.def("save_to", &T::save_to, "Save to file.");
		typedef typename T::eval_type (T::*eval_free)(const double &) const;
		inst.def("eval", eval_free(&T::eval), "Evaluate at time arg2.");
		typedef double (T::*norm_named)() const;
		inst.add_property("norm", norm_named(&T::norm), "Norm.");
		inst.add_property("atoms", &T::atoms, "Number of atoms composing the series.");
		inst.def("swap", &T::swap, "Swap contents with series arg2.");
		// Equality.
		inst.def(boost::python::self == double());
		inst.def(boost::python::self == piranha::mp_rational());
		inst.def(boost::python::self == piranha::mp_integer());
		inst.def(boost::python::self == boost::python::self);
		inst.def(boost::python::self != double());
		inst.def(boost::python::self != piranha::mp_rational());
		inst.def(boost::python::self != piranha::mp_integer());
		inst.def(boost::python::self != boost::python::self);
		// Addition.
		inst.def(boost::python::self += double());
		inst.def(boost::python::self += piranha::mp_rational());
		inst.def(boost::python::self += piranha::mp_integer());
		inst.def(boost::python::self += boost::python::self);
		inst.def(boost::python::self + double());
		inst.def(boost::python::self + piranha::mp_rational());
		inst.def(boost::python::self + piranha::mp_integer());
		inst.def(double() + boost::python::self);
		inst.def(piranha::mp_rational() + boost::python::self);
		inst.def(piranha::mp_integer() + boost::python::self);
		inst.def(boost::python::self + boost::python::self);
		// Subtraction (same as above).
		inst.def(boost::python::self -= double());
		inst.def(boost::python::self -= piranha::mp_rational());
		inst.def(boost::python::self -= piranha::mp_integer());
		inst.def(boost::python::self -= boost::python::self);
		inst.def(boost::python::self - double());
		inst.def(boost::python::self - piranha::mp_rational());
		inst.def(boost::python::self - piranha::mp_integer());
		inst.def(double() - boost::python::self);
		inst.def(piranha::mp_rational() - boost::python::self);
		inst.def(piranha::mp_integer() - boost::python::self);
		inst.def(boost::python::self - boost::python::self);
		// Negation.
		inst.def(-boost::python::self);
		// Multiplication.
		inst.def(boost::python::self *= double());
		inst.def(boost::python::self *= piranha::mp_rational());
		inst.def(boost::python::self *= piranha::mp_integer());
		inst.def(boost::python::self *= boost::python::self);
		inst.def(boost::python::self * double());
		inst.def(boost::python::self * piranha::mp_rational());
		inst.def(boost::python::self * piranha::mp_integer());
		inst.def(double() * boost::python::self);
		inst.def(piranha::mp_rational() * boost::python::self);
		inst.def(piranha::mp_integer() * boost::python::self);
		inst.def(boost::python::self * boost::python::self);
		// Division.
		inst.def(boost::python::self /= double());
		inst.def(boost::python::self /= piranha::mp_rational());
		inst.def(boost::python::self /= piranha::mp_integer());
		inst.def(boost::python::self / double());
		inst.def(boost::python::self / piranha::mp_rational());
		inst.def(boost::python::self / piranha::mp_integer());
		// Exponentiation.
		typedef T (T::*pow_double)(const double &) const;
		typedef T (T::*pow_rational)(const piranha::mp_rational &) const;
		inst.def("__pow__", pow_double(&T::pow));
		inst.def("__pow__", pow_rational(&T::pow));
		inst.def("root", &T::root, "arg2-th root.");
		return std::make_pair(inst, term_inst);
	}

	template <class T>
	inline void series_complex_instantiation(boost::python::class_<std::complex<T> > &instc,
		boost::python::class_<T> &inst)
	{
		// Ctors.
		instc.def(boost::python::init<const std::complex<double> &>());
		instc.def(boost::python::init<const T &>());
		instc.def(boost::python::init<const T &, const T &>());
		inst.def("complex", &T::complex, "Return complex counterpart.");
		instc.def(boost::python::self == T());
		instc.def(boost::python::self != T());
		instc.def(boost::python::self == std::complex<double>());
		instc.def(boost::python::self != std::complex<double>());
		// Addition and subtraction.
		instc.def(boost::python::self += std::complex<double>());
		instc.def(boost::python::self += T());
		instc.def(boost::python::self + std::complex<double>());
		instc.def(std::complex<double>() + boost::python::self);
		instc.def(boost::python::self + T());
		instc.def(T() + boost::python::self);
		instc.def(boost::python::self -= std::complex<double>());
		instc.def(boost::python::self -= T());
		instc.def(boost::python::self - std::complex<double>());
		instc.def(std::complex<double>() - boost::python::self);
		instc.def(boost::python::self - T());
		instc.def(T() - boost::python::self);
		// Multiplication.
		instc.def(boost::python::self *= std::complex<double>());
		instc.def(boost::python::self *= T());
		instc.def(boost::python::self * std::complex<double>());
		instc.def(std::complex<double>() * boost::python::self);
		instc.def(boost::python::self * T());
		instc.def(T() * boost::python::self);
		// Division.
		instc.def(boost::python::self /= std::complex<double>());
		instc.def(boost::python::self / std::complex<double>());
		// Real and imaginary parts assignment and extraction.
		instc.add_property("real", &std::complex<T>::real, &std::complex<T>::set_real, "Get/set real part.");
		instc.add_property("imag", &std::complex<T>::imag, &std::complex<T>::set_imag, "Get/set imaginary part.");
		// Absolute values.
		instc.def("__abs__", &std::complex<T>::abs, "Absolute value.");
		instc.def("abs2", &std::complex<T>::abs2, "Squared absolute value.");
		// Conjugate.
		instc.def("conjugate", &std::complex<T>::conjugate, "Complex conjugate.");
	}

	template <class T>
	inline void series_trigonometric_instantiation(boost::python::class_<T> &inst)
	{
		typedef std::complex<T> (T::*named_ei)() const;
		inst.def("ei", named_ei(&T::ei), "Complex exponential.");
		inst.def("cos", &T::cos, "Cosine.");
		inst.def("sin", &T::sin, "Sine.");
		typedef std::complex<T> (*Ynm_named)(const int &, const int &,
			const T &, const T &);
		typedef std::complex<T> (*Ynm_ei_named)(const int &, const int &,
			const T &, const std::complex<T> &, const std::complex<T> &);
		typedef std::complex<T> (*Ynm_wigner)(const int &, const int &,
			const T &, const T &, const T &, const T &, const T &);
		typedef std::complex<T> (*Ynm_ei_wigner)(const int &, const int &,
			const T &, const std::complex<T> &, const std::complex<T> &, const T &, const T &, const T &);
		const char *Ynm_descr = "Non-normalised spherical harmonic.";
		inst.def("Ynm", Ynm_named(&T::Ynm), Ynm_descr);
		inst.def("Ynm", Ynm_ei_named(&T::Ynm), Ynm_descr);
		inst.def("Ynm", Ynm_wigner(&T::Ynm), Ynm_descr);
		inst.def("Ynm", Ynm_ei_wigner(&T::Ynm), Ynm_descr);
		inst.staticmethod("Ynm");
	}

	template <class T>
	static inline T py_series_partial_psym_1(const T &series, const piranha::psym &p)
	{
		return series.partial(p);
	}

	template <class T>
	static inline T py_series_partial_psym_n(const T &series, const piranha::psym &p, const int &n)
	{
		return series.partial(p,n);
	}

	template <class T>
	inline void series_differential_instantiation(boost::python::class_<T> &inst)
	{
		inst.def("partial", &py_series_partial_psym_1<T>, "Partial derivative with respect to psym arg2.");
		inst.def("partial", &py_series_partial_psym_n<T>, "Partial derivative of natural order arg3 "
			"with respect to psym arg2.");
	}

	template <class T, class Series>
	static inline T py_series_sub_psym_series(const T &series, const piranha::psym &p, const Series &sub)
	{
		return series.template sub<Series>(p, sub);
	}

	template <class T, class Series>
	static inline T py_series_ei_sub_psym_series(const T &series, const piranha::psym &p, const Series &sub)
	{
		return series.template ei_sub<Series>(p, sub);
	}

	template <class T, class Series>
	inline void series_sub_instantiation(boost::python::class_<T> &inst)
	{
		inst.def("sub", &py_series_sub_psym_series<T, Series>, "Substitute psym arg2 with series arg3.");
	}

	template <class T, class Series>
	inline void series_ei_sub_instantiation(boost::python::class_<T> &inst)
	{
		inst.def("ei_sub", &py_series_ei_sub_psym_series<T, Series>, "Substitute complex exponential of psym "
			"arg2 with series arg3 in trigonometric keys.");
	}

	template <class T>
	inline void series_special_functions_instantiation(boost::python::class_<T> &inst)
	{
		typedef T (T::*named_1)(const int &) const;
		typedef T (T::*named_2)(const int &, const int &) const;
		inst.def("besselJ", named_1(&T::besselJ), "Bessel function of the first kind of integer order arg2.");
		inst.def("dbesselJ", named_1(&T::dbesselJ), "Partial derivative of Bessel function of the first kind "
			"of integer order arg2.");
		inst.def("besselJ_div_m", named_2(&T::besselJ_div_m),
			"Bessel function of the first kind of integer order arg2 divided by its argument**arg3.");
		typedef T (T::*named_3)(const int &, const int &, const T &) const;
		const char *Pnm_descr = "Associated Legendre function of integer degree arg2 and order arg3.";
		inst.def("Pnm", named_2(&T::Pnm), Pnm_descr);
		inst.def("Pnm", named_3(&T::Pnm), Pnm_descr);
		inst.def("Pn", &T::Pn, "Legendre polynomial of degree arg2.");
		inst.def("hyperF", &T::template hyperF<piranha::mp_rational>);
	}

#define __celmec_inst(arg) \
	inst.def(#arg,&T::arg,arg##_docstring) \
	.staticmethod(#arg);

	template <class T>
	inline void celmec_instantiation(boost::python::class_<T> &inst)
	{
		const char *r_a_docstring = "Elliptic expansion of r / a in terms of eccentricity arg2 and mean anomaly arg3.";
		const char *a_r_docstring = "Elliptic expansion of a / r in terms of eccentricity arg2 and mean anomaly arg3.";
		const char *cos_f_docstring = "Elliptic expansion of cos(f) in terms of eccentricity arg2 and mean anomaly arg3.";
		const char *sin_f_docstring = "Elliptic expansion of sin(f) in terms of eccentricity arg2 and mean anomaly arg3.";
		const char *cos_E_docstring = "Elliptic expansion of cos(E) in terms of eccentricity arg2 and mean anomaly arg3.";
		const char *sin_E_docstring = "Elliptic expansion of sin(E) in terms of eccentricity arg2 and mean anomaly arg3.";
		const char *EE_docstring = "Elliptic expansion of E (solution of Kepler's equation) "
			"in terms of eccentricity arg2 and mean anomaly arg3.";
		const char *eipE_docstring = "Elliptic expansion of exp(i*p*E) in terms of eccentricity arg2 and mean anomaly arg3, "
			"with p = arg4.";
		__celmec_inst(r_a);
		__celmec_inst(a_r);
		__celmec_inst(cos_f);
		__celmec_inst(sin_f);
		__celmec_inst(cos_E);
		__celmec_inst(sin_E);
		__celmec_inst(EE);
		__celmec_inst(eipE);
	}

#undef __celmec_inst

	template <class T>
	inline void power_series_instantiation(boost::python::class_<T> &inst)
	{
		inst.def("degree", &T::degree, "Degree.");
		inst.def("partial_degree", &T::partial_degree, "Partial degree.");
		inst.def("min_degree", &T::min_degree, "Minimum degree.");
		inst.def("partial_min_degree", &T::partial_min_degree, "Partial minimum degree.");
	}

	template <class T>
	inline void common_polynomial_instantiation(boost::python::class_<T> &inst)
	{
		series_differential_instantiation(inst);
		series_special_functions_instantiation(inst);
		power_series_instantiation(inst);
	}

	template <class T>
	inline void common_poisson_series_instantiation(boost::python::class_<T> &inst, const std::string &)
	{
		series_differential_instantiation(inst);
		series_special_functions_instantiation(inst);
		power_series_instantiation(inst);
		// Expose the polynomial coefficient.
		typedef typename T::term_type::cf_type cf_type;
		//boost::python::class_<cf_type> poly_cf_inst(cf_bindings<cf_type>((name + "_cf").c_str(), ""));
		//power_series_instantiation(poly_cf_inst);
		//poly_cf_inst.def("__iter__", boost::python::iterator<cf_type>());
		//poly_cf_inst.def("__append__", &py_series_append<cf_type,typename cf_type::term_type>);
		//py_series_term<typename cf_type::term_type>(name + "_cf",name);
	}

	template <class T>
	inline void common_fourier_series_instantiation(boost::python::class_<T> &inst)
	{
		series_special_functions_instantiation(inst);
		series_differential_instantiation(inst);
	}
}

#endif
