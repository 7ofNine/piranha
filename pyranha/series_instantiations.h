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
#include <vector>

#include "../src/core/base_classes/named_series_def.h" // For eval_dict.
#include "../src/core/config.h"
#include "../src/core/mp.h"
#include "../src/core/psym.h"
#include "../src/core/type_traits.h"
#include "boost_python_container_conversions.h"
#include "commons.h"

namespace pyranha
{
	// "Fat" iterator that wraps the iterator of named series, allowing retrieval of coefficient and key
	// of each term as a 2-vector of named series. Provides safe and useful iterator capabilities for series in Python.
	// NOT thread-safe.
	template <class NamedSeries>
	struct __PIRANHA_VISIBLE named_iterator {
		typedef typename NamedSeries::const_iterator s_iterator;
		typedef typename s_iterator::iterator_category iterator_category;
		typedef std::vector<NamedSeries> value_type;
		typedef typename s_iterator::difference_type difference_type;
		typedef value_type * pointer;
		typedef value_type & reference;
		named_iterator(const s_iterator &it, const NamedSeries &s):m_it(it),m_s(&s) {}
		named_iterator &operator=(const named_iterator &n)
		{
			if (*this != &n) {
				m_it = n.m_it;
				m_s = n.m_s;
			}
			return *this;
		}
		value_type &operator*() const
		{
			if (m_value.size() != 2) {
				m_value.resize(2);
			}
			m_value[0] = m_s->series_from_cf(m_it->m_cf);
			m_value[1] = m_s->series_from_key(m_it->m_key);
			return m_value;
		}
		bool operator==(const named_iterator &n_it) const
		{
			return (m_it == n_it.m_it && m_s == n_it.m_s);
		}
		named_iterator &operator++()
		{
			++m_it;
			return *this;
		}
		named_iterator operator++(int)
		{
			named_iterator retval(m_it,*m_s);
			++m_it;
			return retval;
		}
		s_iterator 		m_it;
		const NamedSeries	*m_s;
		static value_type	m_value;
	};

	template <class NamedSeries>
	inline named_iterator<NamedSeries> series_begin(const NamedSeries &s)
	{
		return named_iterator<NamedSeries>(s.begin(),s);
	}

	template <class NamedSeries>
	inline named_iterator<NamedSeries> series_end(const NamedSeries &s)
	{
		return named_iterator<NamedSeries>(s.end(),s);
	}

	// Initializer for named_iterator's static member.
	template <class NamedSeries>
	std::vector<NamedSeries> named_iterator<NamedSeries>::m_value;

	/// Basic series instantiation.
	template <class T>
	inline boost::python::class_<T> series_basic_instantiation(const std::string &name, const std::string &description)
	{
		// Expose vectors of series and vectors of vectors of series as Python tuples.
		to_tuple_mapping<std::vector<T> >();
		to_tuple_mapping<std::vector<std::vector<T> > >();
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
		//inst.def("__iter__", boost::python::iterator<T>());
		inst.def("__iter__", boost::python::range(&series_begin<T>, &series_end<T>));
		inst.def("__len__", &T::length);
		inst.def("__impl_repr__", &py_print_to_string<T>);
		// Pyranha-specific special methods.
		inst.def("_latex_", &py_print_to_string_tex<T>, "Latex representation.");
		inst.def("dump", &py_print_to_string_plain<T>, "Return a string of the series in plain format.");
		inst.add_property("arguments", &py_series_arguments<T>, "Series arguments.");
		inst.def("__split__", &T::split, "Split series.");
		inst.def("__psi__", &T::psi, "Power series iterations.");
		inst.def("save_to", &T::save_to, "Save to file.");
		typedef typename T::eval_type (T::*eval_double)(const double &) const;
		typedef typename T::eval_type (T::*eval_dic)(const piranha::eval_dict &) const;
		inst.def("__eval__", eval_double(&T::eval), "Evaluate at time arg2.");
		inst.def("__eval__", eval_dic(&T::eval), "Evaluate using dictionary arg2.");
		inst.add_property("norm", &T::norm, "Norm.");
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
		return inst;
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
		inst.def("ei", &T::ei, "Complex exponential.");
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
	static inline T py_series_partial_1(const T &series, const std::string &s)
	{
		return series.partial(s);
	}

	template <class T>
	static inline T py_series_partial_n(const T &series, const std::string &s, const int &n)
	{
		return series.partial(s,n);
	}

	template <class T>
	inline void series_differential_instantiation(boost::python::class_<T> &inst)
	{
		inst.def("partial", &py_series_partial_1<T>, "Partial derivative with respect to psym named arg2.");
		inst.def("partial", &py_series_partial_n<T>, "Partial derivative of natural order arg3 "
			"with respect to psym named arg2.");
	}

	template <class T>
	inline void series_integral_instantiation(boost::python::class_<T> &inst)
	{
		inst.def("integrate", &T::integrate, "Indefinite integration with respect to psym named arg2.");
	}

	template <class T, class Series>
	static inline T py_series_sub_string_series(const T &series, const std::string &s, const Series &sub)
	{
		return series.template sub<Series>(s, sub);
	}

	template <class T, class Series>
	static inline T py_series_ei_sub_string_series(const T &series, const std::string &s, const Series &sub)
	{
		return series.template ei_sub<Series>(s, sub);
	}

	template <class T, class Series>
	inline void series_sub_instantiation(boost::python::class_<T> &inst)
	{
		inst.def("sub", &py_series_sub_string_series<T, Series>, "Substitute psym named arg2 with series arg3.");
	}

	template <class T, class Series>
	inline void series_ei_sub_instantiation(boost::python::class_<T> &inst)
	{
		inst.def("ei_sub", &py_series_ei_sub_string_series<T, Series>, "Substitute complex exponential of psym "
			"named arg2 with series arg3 in trigonometric keys.");
	}

	template <class T>
	inline void series_special_functions_instantiation(boost::python::class_<T> &inst)
	{
		inst.def("besselJ", &T::besselJ, "Bessel function of the first kind of integer order arg2.");
		inst.def("dbesselJ", &T::dbesselJ, "Partial derivative of Bessel function of the first kind "
			"of integer order arg2.");
		inst.def("besselJ_div_m", &T::besselJ_div_m,
			"Bessel function of the first kind of integer order arg2 divided by its argument**arg3.");
		typedef T (T::*named_2)(const int &, const int &) const;
		typedef T (T::*named_3)(const int &, const int &, const T &) const;
		const char *Pnm_descr = "Associated Legendre function of integer degree arg2 and order arg3.";
		inst.def("legendrePnm", named_2(&T::legendrePnm), Pnm_descr);
		inst.def("legendrePnm", named_3(&T::legendrePnm), Pnm_descr);
		inst.def("legendrePn", &T::legendrePn, "Legendre polynomial of degree arg2.");
		typedef T (T::*hyper_1)(const std::vector<piranha::mp_rational> &, const std::vector<piranha::mp_rational> &) const;
		typedef T (T::*hyper_2)(const std::vector<piranha::mp_rational> &, const std::vector<piranha::mp_rational> &, const int &) const;
		inst.def("hyperF", hyper_1(&T::hyperF));
		inst.def("hyperF", hyper_2(&T::hyperF));
		typedef T (T::*dhyper_1)(const int &, const std::vector<piranha::mp_rational> &, const std::vector<piranha::mp_rational> &) const;
		typedef T (T::*dhyper_2)(const int &, const std::vector<piranha::mp_rational> &, const std::vector<piranha::mp_rational> &, const int &) const;
		inst.def("dhyperF", dhyper_1(&T::dhyperF));
		inst.def("dhyperF", dhyper_2(&T::dhyperF));
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

	template <class PowerSeries>
	typename PowerSeries::degree_type p_degree_str(const PowerSeries &p, const std::string &name)
	{
		return p.partial_degree(std::vector<std::string>(1,name));
	}

	template <class PowerSeries>
	typename PowerSeries::degree_type p_order_str(const PowerSeries &p, const std::string &name)
	{
		return p.partial_order(std::vector<std::string>(1,name));
	}

	template <class T>
	inline void power_series_instantiation(boost::python::class_<T> &inst)
	{
		inst.def("degree", &T::degree, "(Partial) degree.");
		inst.def("degree", &T::partial_degree);
		inst.def("degree", &p_degree_str<T>);
		inst.def("order", &T::order, "(Partial) order.");
		inst.def("order", &T::partial_order);
		inst.def("order", &p_order_str<T>);
	}

	template <class HarmonicSeries>
	typename HarmonicSeries::h_degree_type p_h_degree_str(const HarmonicSeries &p, const std::string &name)
	{
		return p.partial_h_degree(std::vector<std::string>(1,name));
	}

	template <class HarmonicSeries>
	typename HarmonicSeries::h_degree_type p_h_order_str(const HarmonicSeries &p, const std::string &name)
	{
		return p.partial_h_order(std::vector<std::string>(1,name));
	}

	template <class T>
	inline void harmonic_series_instantiation(boost::python::class_<T> &inst)
	{
		inst.def("h_degree", &T::h_degree, "(Partial) harmonic degree.");
		inst.def("h_degree", &T::partial_h_degree);
		inst.def("h_degree", &p_h_degree_str<T>);
		inst.def("h_order", &T::h_order, "(Partial) harmonic order.");
		inst.def("h_order", &T::partial_h_order);
		inst.def("h_order", &p_h_order_str<T>);
		inst.def("is_sine", &T::is_sine, "Return true if series is made only of sine terms.");
		inst.def("is_cosine", &T::is_cosine, "Return true if series is made only of cosine terms.");
	}

	template <class T>
	inline void common_polynomial_instantiation(boost::python::class_<T> &inst)
	{
		series_differential_instantiation(inst);
		series_integral_instantiation(inst);
		series_special_functions_instantiation(inst);
		power_series_instantiation(inst);
	}

	template <class T>
	inline void common_poisson_series_instantiation(boost::python::class_<T> &inst, const std::string &)
	{
		series_differential_instantiation(inst);
		series_special_functions_instantiation(inst);
		power_series_instantiation(inst);
		harmonic_series_instantiation(inst);
		series_integral_instantiation(inst);
	}

	template <class T>
	inline void common_fourier_series_instantiation(boost::python::class_<T> &inst)
	{
		series_special_functions_instantiation(inst);
		series_differential_instantiation(inst);
		harmonic_series_instantiation(inst);
		series_integral_instantiation(inst);
	}
}

#endif
