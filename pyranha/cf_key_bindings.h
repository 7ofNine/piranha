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
#include <vector>

#include "../src/core/exceptions.h"
#include "../src/core/ntuple.h"
#include "../src/core/numerical_coefficients/double_cf.h"
#include "../src/core/numerical_coefficients/mpq_cf.h"
#include "../src/core/numerical_coefficients/mpz_cf.h"
#include "../src/core/poisson_series_common/trig_array.h"
#include "../src/core/polynomial_common/expo_array.h"
#include "../src/core/psym.h"
#include "commons.h"

namespace pyranha
{
	// Secure wrappers to access methods of cfkeys that require an args tuple as input.
	template <class CfKey, class ArgsTuple>
	static inline double py_cfkey_norm(const CfKey &cfkey, const ArgsTuple &a)
	{
		if (!cfkey.is_insertable(a)) {
			throw piranha::unsuitable("Arguments tuple is not compatible with the calculation of the norm "
				"of this coefficient/key.");
		}
		return cfkey.norm_(a);
	}

	template <class CfKey, class ArgsTuple>
	static inline typename CfKey::eval_type py_cfkey_eval(const CfKey &cfkey, const double &t, const ArgsTuple &a)
	{
		if (!cfkey.is_insertable(a)) {
			throw piranha::unsuitable("Arguments tuple is not compatible with the calculation "
				"of the time evaluation of this coefficient/key.");
		}
		return cfkey.eval_(t,a);
	}

	template <class Cf, class ArgsTuple>
	static inline Cf py_cf_real(const std::complex<Cf> &cf, const ArgsTuple &a)
	{
		if (!cf.is_insertable(a)) {
			throw piranha::unsuitable("Arguments tuple is not compatible with the extraction "
				"of the real part from this coefficient.");
		}
		return cf.real_(a);
	}

	template <class Cf, class ArgsTuple>
	static inline Cf py_cf_imag(const std::complex<Cf> &cf, const ArgsTuple &a)
	{
		if (!cf.is_insertable(a)) {
			throw piranha::unsuitable("Arguments tuple is not compatible with the extraction "
				"of the imaginary part from this coefficient.");
		}
		return cf.imag_(a);
	}

	// Bindings common to coefficients and keys that require args tuple as input. Instantiations for
	// args tuples of different sizes.
	template <int N, class CfKey, int M>
	struct mp_cfkey_bindings {
		static void run(boost::python::class_<CfKey> &cfkey_inst) {
			typedef typename piranha::ntuple<std::vector<piranha::psym_p>,N>::type args_tuple_type;
			cfkey_inst.def("norm", &py_cfkey_norm<CfKey,args_tuple_type>, "Norm.");
			cfkey_inst.def("eval", &py_cfkey_eval<CfKey,args_tuple_type>, "Time evaluation.");
			mp_cfkey_bindings<N - 1,CfKey,M>::run(cfkey_inst);
		}
	};

	template <int N, class CfKey>
	struct mp_cfkey_bindings<N,CfKey,N> {
		static void run(boost::python::class_<CfKey> &) {}
	};

	// Bindings common to coefficients and keys.
	template <class CfKey, int M = 0>
	struct cfkey_bindings {
		static void run(boost::python::class_<CfKey> &cfkey_inst)
		{
			cfkey_inst.def(boost::python::init<const CfKey &>());
			cfkey_inst.def("__copy__",&py_copy<CfKey>);
			cfkey_inst.def("atoms", &CfKey::atoms, "Number of atoms.");
			mp_cfkey_bindings<__PIRANHA_MAX_ECHELON_LEVEL + 1, CfKey, M>::run(cfkey_inst);
		}
	};

	template <int N, class Cf>
	struct mp_complex_cf_bindings {
		static void run(boost::python::class_<std::complex<Cf> > &cf_inst) {
			typedef typename piranha::ntuple<std::vector<piranha::psym_p>,N>::type args_tuple_type;
			cf_inst.def("real", &py_cf_real<Cf,args_tuple_type>, "Real part.");
			cf_inst.def("imag", &py_cf_imag<Cf,args_tuple_type>, "Imaginary part.");
			mp_complex_cf_bindings<N - 1,Cf>::run(cf_inst);
		}
	};

	template <class Cf>
	struct mp_complex_cf_bindings<0,Cf> {
		static inline void run(boost::python::class_<std::complex<Cf> > &) {}
	};

	template <class Cf>
	static inline void complex_cf_bindings(boost::python::class_<std::complex<Cf> > &cf_inst)
	{
		mp_complex_cf_bindings<__PIRANHA_MAX_ECHELON_LEVEL + 1, Cf>::run(cf_inst);
	}

	inline void numerical_cfs_bindings()
	{
		using namespace boost::python;
		using namespace piranha;
		using namespace std;
		class_<double_cf> dcf("double_cf", "Double precision coefficient.");
		class_<mpz_cf> zcf("mpz_cf", "Arbitrary precision integer coefficient.");
		class_<mpq_cf> qcf("mpq_cf", "Arbitrary precision rational coefficient.");
		cfkey_bindings<double_cf>::run(dcf);
		cfkey_bindings<mpz_cf>::run(zcf);
		cfkey_bindings<mpq_cf>::run(qcf);
		// Now the complex ones.
		class_<complex<double_cf> > cdcf("double_cfc", "Complex double precision coefficient.");
		class_<complex<mpz_cf> > czcf("mpz_cfc", "Complex arbitrary precision integer coefficient.");
		class_<complex<mpq_cf> > cqcf("mpq_cfc", "Complex arbitrary precision rational coefficient.");
		cfkey_bindings<complex<double_cf> >::run(cdcf);
		cfkey_bindings<complex<mpz_cf> >::run(czcf);
		cfkey_bindings<complex<mpq_cf> >::run(cqcf);
		// Complex-specific stuff.
		complex_cf_bindings<double_cf>(cdcf);
		complex_cf_bindings<mpz_cf>(czcf);
		complex_cf_bindings<mpq_cf>(cqcf);
	}

	template <class TrigArray, class ArgsTuple>
	static inline double py_trigarray_freq(const TrigArray &t, const ArgsTuple &a)
	{
		if (!t.is_insertable(a)) {
			throw piranha::unsuitable("Arguments tuple is not compatible with the calculation of the frequency "
				"of this trigonometric array.");
		}
		return t.freq(a);
	}

	template <class TrigArray, class ArgsTuple>
	static inline double py_trigarray_phase(const TrigArray &t, const ArgsTuple &a)
	{
		if (!t.is_insertable(a)) {
			throw piranha::unsuitable("Arguments tuple is not compatible with the calculation of the phase "
				"of this trigonometric array.");
		}
		return t.phase(a);
	}

	template <int N, class TrigArray, int M>
	struct mp_ta_bindings {
		static void run(boost::python::class_<TrigArray> &ta_inst) {
			typedef typename piranha::ntuple<std::vector<piranha::psym_p>,N>::type args_tuple_type;
			ta_inst.def("freq", &py_trigarray_freq<TrigArray,args_tuple_type>, "Frequency.");
			ta_inst.def("phase", &py_trigarray_phase<TrigArray,args_tuple_type>, "Phase.");
			mp_ta_bindings<N - 1,TrigArray, M>::run(ta_inst);
		}
	};

	template <int N, class TrigArray>
	struct mp_ta_bindings<N,TrigArray,N> {
		static void run(boost::python::class_<TrigArray> &) {}
	};

	template <class TrigArray, int M = 0>
	struct trig_array_bindings {
		static void run(boost::python::class_<TrigArray> &inst)
		{
			inst.add_property("flavour", &TrigArray::get_flavour);
			mp_ta_bindings<__PIRANHA_MAX_ECHELON_LEVEL + 1,TrigArray,M>::run(inst);
		}
	};

	template <class ExpoArray>
	static inline void expo_array_bindings(boost::python::class_<ExpoArray> &inst)
	{
		inst.add_property("degree", &ExpoArray::degree);
	}

	inline void keys_bindings()
	{
		// TODO: re-establish int_array_bindings (for, e.g., __getitem__).
		using namespace boost::python;
		using namespace piranha;
		class_<trig_array<16,0> > ta_16_0("trig_array_16_0");
		class_<trig_array<16,1> > ta_16_1("trig_array_16_1");
		class_<expo_array<16,0> > ea_16_0("expo_array_16_0");
		// Common stuff.
		cfkey_bindings<trig_array<16,0> >::run(ta_16_0);
		cfkey_bindings<trig_array<16,1>,1>::run(ta_16_1);
		cfkey_bindings<expo_array<16,0> >::run(ea_16_0);
		// Trig-array specifics.
		trig_array_bindings<trig_array<16,0> >::run(ta_16_0);
		trig_array_bindings<trig_array<16,1>,1>::run(ta_16_1);
		// Expo-array specifics.
		expo_array_bindings(ea_16_0);
	}
}

#endif
