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

#include "../src/core/numerical_coefficients/double_cf.h"
#include "../src/core/numerical_coefficients/mpq_cf.h"
#include "../src/core/numerical_coefficients/mpz_cf.h"
#include "../src/core/poisson_series_common/trig_array.h"
#include "../src/core/shared_args.h"

namespace pyranha
{
	// TODO: Do _not_ provide a free or py_interface in classes, use wrappers here.
	// In the wrappers error checks can be performed.
	template <class Cf>
	inline boost::python::class_<Cf> cf_bindings(const std::string &name, const std::string &description)
	{
		boost::python::class_<Cf> cf_inst(name.c_str(), description.c_str());
		typedef double(Cf::*norm_free)() const;
		cf_inst.def("norm", norm_free(&Cf::norm), "Norm.");
		cf_inst.def("atoms", &Cf::atoms, "Number of atoms.");
		return cf_inst;
	}

	template <class Cf>
	inline boost::python::class_<std::complex<Cf> > complex_cf_bindings(const std::string &name,
			const std::string &description)
	{
		boost::python::class_<std::complex<Cf> > cf_inst(cf_bindings<std::complex<Cf> >(name, description));
		typedef Cf(std::complex<Cf>::*comp_free)() const;
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

	// TODO: this will be moved in a common header later.
	template <class Vector>
	inline typename Vector::value_type py_vector_getitem(const Vector &v, const piranha::max_fast_int &n_)
	{
		int n = n_;
		const size_t size = v.size();
		if (n_ < 0) {
			n = n_ + size;
		}
		if (n < 0 || static_cast<size_t>(n) >= size) {
			throw piranha::unsuitable("Index of integer array key is out of range");
		}
		return v[n];
	}

	template <class CfKey>
	inline double py_cfkey_norm(const CfKey &cfkey)
	{
		return cfkey.norm(piranha::shared_args::get());
	}

	template <class IntArrayKey>
	inline boost::python::class_<IntArrayKey> int_array_key_bindings(const std::string &name, const std::string &descr)
	{
		boost::python::class_<IntArrayKey> retval(name.c_str(), descr.c_str());
		retval.def("__getitem__", &py_vector_getitem<IntArrayKey>);
		retval.def("norm", py_cfkey_norm<IntArrayKey>);
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
		return t.freq(piranha::shared_args::get());
	}

	template <class TrigArray>
	inline double py_trigarray_phase(const TrigArray &t)
	{
		return t.phase(piranha::shared_args::get());
	}

	template <class TrigArray>
	inline void trig_array_key_bindings(boost::python::class_<TrigArray> &inst)
	{
		inst.add_property("flavour", &py_trigarray_flavour<TrigArray>);
		inst.def("freq", &py_trigarray_freq<TrigArray>);
		inst.def("phase", &py_trigarray_phase<TrigArray>);
	}

	inline void keys_bindings()
	{
		boost::python::class_<piranha::trig_array<16, 0> >
		ta_16_0(int_array_key_bindings<piranha::trig_array<16, 0> >("trig_array_16_0", ""));
		trig_array_key_bindings(ta_16_0);

		boost::python::class_<piranha::trig_array<16, 1> >
		ta_16_1(int_array_key_bindings<piranha::trig_array<16, 1> >("trig_array_16_1", ""));
		trig_array_key_bindings(ta_16_1);
	}
}

#endif
