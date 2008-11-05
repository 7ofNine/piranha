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

#ifndef PYRANHA_CF_KEY_BINDINGS_H
#define PYRANHA_CF_KEY_BINDINGS_H

#include <complex>
#include <boost/python/class.hpp>
#include <string>

#include "../src/core/exceptions.h"
#include "../src/core/numerical_coefficients/double_cf.h"
#include "../src/core/numerical_coefficients/mpq_cf.h"
#include "../src/core/numerical_coefficients/mpz_cf.h"
#include "../src/core/poisson_series_common/trig_array.h"
#include "../src/core/polynomial_common/expo_array.h"
#include "../src/core/shared_args.h"
#include "commons.h"

namespace pyranha
{
	template <class CfKey>
	inline double py_cfkey_norm(const CfKey &cfkey)
	{
		if (!cfkey.is_insertable(piranha::shared_args::get())) {
			throw piranha::unsuitable("Current shared arguments are not compatible with the calculation of the norm "
				"of this coefficient/key.");
		}
		return cfkey.norm(piranha::shared_args::get());
	}

	template <class CfKey>
	inline typename CfKey::eval_type py_cfkey_eval(const CfKey &cfkey, const double &t)
	{
		if (!cfkey.is_insertable(piranha::shared_args::get())) {
			throw piranha::unsuitable("Current shared arguments are not compatible with the calculation "
				"of the time evaluation of this coefficient/key.");
		}
		return cfkey.eval(t,piranha::shared_args::get());
	}

	// TODO: Do _not_ provide a free or py_interface in classes, use wrappers here.
	// In the wrappers error checks can be performed.
	template <class Cf>
	inline boost::python::class_<Cf> cf_bindings(const std::string &name, const std::string &description)
	{
		boost::python::class_<Cf> cf_inst(name.c_str(), description.c_str());
		cf_inst.def(boost::python::init<const Cf &>());
		cf_inst.def("__abs__", &py_cfkey_norm<Cf>);
		cf_inst.def("__copy__",&py_copy<Cf>);
		cf_inst.add_property("norm", &py_cfkey_norm<Cf>, "Norm.");
		cf_inst.def("eval", &py_cfkey_eval<Cf>, "Time evaluation.");
		cf_inst.add_property("atoms", &Cf::atoms, "Number of atoms.");
		return cf_inst;
	}

	template <class Cf>
	inline boost::python::class_<std::complex<Cf> > complex_cf_bindings(const std::string &name,
			const std::string &description)
	{
		boost::python::class_<std::complex<Cf> > cf_inst(cf_bindings<std::complex<Cf> >(name, description));
		typedef Cf(std::complex<Cf>::*comp_free)() const;
		// TODO: ditch this free interface, use wrappers.
		cf_inst.def("real", comp_free(&std::complex<Cf>::real), "Real part.");
		cf_inst.def("imag", comp_free(&std::complex<Cf>::imag), "Imaginary part.");
		return cf_inst;
	}

	inline void numerical_cfs_bindings()
	{
		cf_bindings<piranha::double_cf>("double_cf", "Double precision coefficient.");
		complex_cf_bindings<piranha::double_cf>("double_cfc", "Complex double precision coefficient.");
		cf_bindings<piranha::mpq_cf>("mpq_cf", "Arbitrary precision rational coefficient.");
		complex_cf_bindings<piranha::mpq_cf>("mpq_cfc", "Complex arbitrary precision rational coefficient.");
		cf_bindings<piranha::mpz_cf>("mpz_cf", "Arbitrary precision integer coefficient.");
		complex_cf_bindings<piranha::mpz_cf>("mpz_cfc", "Complex arbitrary precision integer coefficient.");
	}

	template <class IntArrayKey>
	inline boost::python::class_<IntArrayKey> int_array_key_bindings(const std::string &name, const std::string &descr)
	{
		boost::python::class_<IntArrayKey> retval(name.c_str(), descr.c_str());
		retval.def("__getitem__", &py_vector_getitem<IntArrayKey>);
		retval.add_property("norm", py_cfkey_norm<IntArrayKey>, "Norm.");
		retval.def("eval", &py_cfkey_eval<IntArrayKey>, "Time evaluation.");
		return retval;
	}

	template <class TrigArray>
	inline bool py_trigarray_flavour(const TrigArray &t)
	{
		return t.flavour();
	}

	template <class TrigArray>
	inline double py_trigarray_freq(const TrigArray &t)
	{
		if (!t.is_insertable(piranha::shared_args::get())) {
			throw piranha::unsuitable("Current shared arguments are not compatible with the calculation of the frequency "
				"of this trigonometric array.");
		}
		return t.freq(piranha::shared_args::get());
	}

	template <class TrigArray>
	inline double py_trigarray_phase(const TrigArray &t)
	{
		if (!t.is_insertable(piranha::shared_args::get())) {
			throw piranha::unsuitable("Current shared arguments are not compatible with the calculation of the phase "
				"of this trigonometric array.");
		}
		return t.phase(piranha::shared_args::get());
	}

	template <class TrigArray>
	inline void trig_array_key_bindings(boost::python::class_<TrigArray> &inst)
	{
		inst.add_property("flavour", &py_trigarray_flavour<TrigArray>);
		inst.add_property("freq", &py_trigarray_freq<TrigArray>);
		inst.add_property("phase", &py_trigarray_phase<TrigArray>);
	}

	template <class ExpoArray>
	inline void expo_array_key_bindings(boost::python::class_<ExpoArray> &inst)
	{
		inst.add_property("degree", &ExpoArray::degree);
		inst.add_property("min_degree", &ExpoArray::min_degree);
	}

	inline void keys_bindings()
	{
		boost::python::class_<piranha::trig_array<16, 0> >
		ta_16_0(int_array_key_bindings<piranha::trig_array<16, 0> >("trig_array_16_0", ""));
		trig_array_key_bindings(ta_16_0);

		boost::python::class_<piranha::trig_array<16, 1> >
		ta_16_1(int_array_key_bindings<piranha::trig_array<16, 1> >("trig_array_16_1", ""));
		trig_array_key_bindings(ta_16_1);

		boost::python::class_<piranha::expo_array<16, 0> >
		ea_16_0(int_array_key_bindings<piranha::expo_array<16, 0> >("expo_array_16_0", ""));
		expo_array_key_bindings(ea_16_0);
	}
}

#endif
