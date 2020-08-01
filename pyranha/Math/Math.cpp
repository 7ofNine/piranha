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

#include <boost/math/special_functions/gamma.hpp>
#include <boost/python/def.hpp>
#include <boost/python/module.hpp>
#include <boost/python/docstring_options.hpp>
#include <complex>

#include "../../src/core/math.h"
#include "../../src/core/mp.h"
#include "../exceptions.h"

using namespace pyranha;
using namespace piranha;
using namespace boost::python;

BOOST_PYTHON_MODULE(_Math)
{
    docstring_options docOptions(true, false, false);
    translate_exceptions();
    typedef std::complex<double> (*jn_complex)(const int &, const std::complex<double> &);
    typedef double (*jn_real)(const int &, const double &);
    def("besselJ", jn_real(&besselJ),
        "Bessel function of the first kind of integer order and real double-precision argument.");
    def("dbesselJ", jn_real(&dbesselJ),
        "First derivative of the Bessel function of the first kind of integer order and real double-precision argument.");
    def("besselI", &besselI, "Modified Bessel function of the first kind of integer order.");
    def("legendrePnm", &legendrePnm, "Associated Legendre function.");
    def("legendrePn", &legendrePn, "Legendre polynomial.");
    typedef std::complex<double> (*Ynm_plain)(const int &, const int &, const double &, const double &);
    typedef std::complex<double> (*Ynm_ei_plain)(const int &, const int &, const double &,
        const std::complex<double> &, const std::complex<double> &);
    typedef std::complex<double> (*Ynm_rot)(const int &, const int &, const double &, const double &,
        const double &, const double &, const double &);
    typedef std::complex<double> (*Ynm_ei_rot)(const int &, const int &, const double &,
        const std::complex<double> &, const std::complex<double> &, const double &, const double &, const double &);
    def("Ynm", Ynm_plain(&Ynm), "Spherical harmonic (not normalised).");
    def("Ynm", Ynm_ei_plain(&Ynm), "Spherical harmonic (not normalised).");
    def("Ynm", Ynm_rot(&Ynm), "Rotated spherical harmonic (not normalised).");
    def("Ynm", Ynm_ei_rot(&Ynm), "Rotated spherical harmonic (not normalised).");
    // Factorial.
    typedef double (*factorial_double)(const int &);
    typedef mp_integer (*factorial_mp)(const mp_integer &);
    def("__factorial", factorial_double(&factorial), "Factorial (double-precision).");
    def("__factorial", factorial_mp(&factorial), "Factorial (arbitrary precision).");
    // Double factorial.
    typedef mp_integer (*double_factorial_mp)(const mp_integer &);
    typedef double (*double_factorial_double)(const int &);
    def("__double_factorial", double_factorial_mp(&double_factorial), "Double factorial of non-negative integer argument.");
    def("__double_factorial", double_factorial_double(&double_factorial), "Double factorial of non-negative integer argument.");
    // Rising factorial.
    typedef double (*r_factorial_double)(const double &, const int &);
    typedef std::complex<double> (*r_factorial_complex_double)(const std::complex<double> &, const int &);
    typedef mp_integer (*r_factorial_mp_int)(const mp_integer &, const int &);
    typedef mp_rational (*r_factorial_mp_rat)(const mp_rational &, const int &);
    def("__r_factorial", r_factorial_complex_double(&r_factorial), "Rising factorial.");
    def("__r_factorial", r_factorial_double(&r_factorial), "Rising factorial.");
    def("__r_factorial", r_factorial_mp_int(&r_factorial), "Rising factorial.");
    def("__r_factorial", r_factorial_mp_rat(&r_factorial), "Rising factorial.");
    // Falling factorial.
    typedef double (*f_factorial_double)(const double &, const int &);
    typedef std::complex<double> (*f_factorial_complex_double)(const std::complex<double> &, const int &);
    typedef mp_integer (*f_factorial_mp_int)(const mp_integer &, const int &);
    typedef mp_rational (*f_factorial_mp_rat)(const mp_rational &, const int &);
    def("__f_factorial", f_factorial_complex_double(&f_factorial), "Falling factorial.");
    def("__f_factorial", f_factorial_double(&f_factorial), "Falling factorial.");
    def("__f_factorial", f_factorial_mp_int(&f_factorial), "Falling factorial.");
    def("__f_factorial", f_factorial_mp_rat(&f_factorial), "Falling factorial.");
    // Choose function.
    typedef double (*choose_double)(const int &, const int &);
    typedef double (*choose_double_double)(const double &, const int &);
    typedef std::complex<double> (*c_choose_double)(const std::complex<int> &, const int &);
    typedef std::complex<double> (*c_choose_double_double)(const std::complex<double> &, const int &);
    typedef mp_integer (*choose_z)(const mp_integer &, const int &);
    typedef mp_rational (*choose_q)(const mp_rational &, const int &);
    def("__choose", c_choose_double_double(&choose), "Binomial coefficient (complex double-precision).");
    def("__choose", c_choose_double(&choose), "Binomial coefficient (complex double-precision).");
    def("__choose", choose_double_double(&choose), "Binomial coefficient (double-precision).");
    def("__choose", choose_double(&choose), "Binomial coefficient (double-precision).");
    def("__choose", choose_z(&choose), "Binomial coefficient (multiprecision integer).");
    def("__choose", choose_q(&choose), "Binomial coefficient (multiprecision rational).");
    typedef double (*double_gamma)(double);
    def("gamma", double_gamma(&boost::math::tgamma<double>), "Gamma function.");
    typedef int (*cs_int)(const int &);
    typedef int (*cs_z)(const mp_integer &);
    def("cs_phase", cs_int(&cs_phase), "Condon-Shortley phase = (-1)**arg1.");
    def("cs_phase", cs_z(&cs_phase), "Condon-Shortley phase = (-1)**arg1.");
}
