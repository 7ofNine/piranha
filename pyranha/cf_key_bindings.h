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

namespace pyranha
{
	template <class Cf>
	inline boost::python::class_<Cf> numerical_cf_bindings(const std::string &name, const std::string &description) {
		boost::python::class_<Cf> cf_inst(name.c_str(),description.c_str());
		typedef double (Cf::*norm_free)() const;
		cf_inst.def("norm",norm_free(&Cf::norm),"Norm.");
		cf_inst.def("atoms",&Cf::atoms,"Number of atoms.");
		return cf_inst;
	}

	template <class Cf>
	inline boost::python::class_<std::complex<Cf> > complex_numerical_cf_bindings(const std::string &name,
		const std::string &description) {
		boost::python::class_<std::complex<Cf> > cf_inst(numerical_cf_bindings<std::complex<Cf> >(name,description));
		typedef Cf (std::complex<Cf>::*comp_free)() const;
		cf_inst.def("real",comp_free(&std::complex<Cf>::real),"Real part.");
		cf_inst.def("imag",comp_free(&std::complex<Cf>::imag),"Imaginary part.");
		return cf_inst;
	}

	inline void numerical_cfs_bindings() {
		numerical_cf_bindings<piranha::double_cf>("double_cf","Double precision coefficient.");
		complex_numerical_cf_bindings<piranha::double_cf>("double_cfc","Complex double precision coefficient.");
		numerical_cf_bindings<piranha::mpq_cf>("mpq_cf","Arbitrary precision rational coefficient.");
		complex_numerical_cf_bindings<piranha::mpq_cf>("mpq_cfc","Complex arbitrary precision rational coefficient.");
		numerical_cf_bindings<piranha::mpz_cf>("mpz_cf","Arbitrary precision integer coefficient.");
		complex_numerical_cf_bindings<piranha::mpz_cf>("mpz_cfc","Complex arbitrary precision integer coefficient.");
	}
}

#endif
