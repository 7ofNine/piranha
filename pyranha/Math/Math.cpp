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


#include "../../src/core/math.h"
#include "../../src/core/mp.h"
#include "../exceptions.h"

#include "pybind11/pybind11.h"
#include "pybind11/complex.h"
#include "pybind11/stl.h"
#include <complex>


PYBIND11_MODULE(_Math, m)
{
    //docstring_options docOptions(true, false, false);
    pyranha::translate_exceptions();
    typedef std::complex<double> (*jn_complex)(const int &, const std::complex<double> &);
    typedef double (*jn_real)(const int &, const double &);
    m.def("besselJ", jn_real(&piranha::besselJ),
            "Bessel function of the first kind of integer order and real double-precision argument.");
    m.def("dbesselJ", jn_real(&piranha::dbesselJ),
            "First derivative of the Bessel function of the first kind of integer order and real double-precision argument.");
    m.def("besselI", &piranha::besselI, "Modified Bessel function of the first kind of integer order.");
    m.def("legendrePnm", &piranha::legendrePnm, "Associated Legendre function.");
    m.def("legendrePn", &piranha::legendrePn, "Legendre polynomial.");

    typedef std::complex<double> (*Ynm_plain)(const int &, const int &, const double &, const double &);
    typedef std::complex<double> (*Ynm_ei_plain)(const int &, const int &, const double &,
        const std::complex<double> &, const std::complex<double> &);
    typedef std::complex<double> (*Ynm_rot)(const int &, const int &, const double &, const double &,
        const double &, const double &, const double &);
    typedef std::complex<double> (*Ynm_ei_rot)(const int &, const int &, const double &,
        const std::complex<double> &, const std::complex<double> &, const double &, const double &, const double &);
    m.def("Ynm", Ynm_plain(&piranha::Ynm), "Spherical harmonic (not normalised).");
    m.def("Ynm", Ynm_ei_plain(&piranha::Ynm), "Spherical harmonic (not normalised).");
    m.def("Ynm", Ynm_rot(&piranha::Ynm), "Rotated spherical harmonic (not normalised).");
    m.def("Ynm", Ynm_ei_rot(&piranha::Ynm), "Rotated spherical harmonic (not normalised).");
    
    // Factorial.
    typedef double (*factorial_double)(const int &);
    typedef piranha::mp_integer (*factorial_mp)(const piranha::mp_integer &);
    m.def("__factorial", factorial_double(&piranha::factorial), "Factorial (double-precision).");
    m.def("__factorial", factorial_mp(&piranha::factorial), "Factorial (arbitrary precision).");
    
    // Double factorial.
    typedef piranha::mp_integer (*double_factorial_mp)(const piranha::mp_integer &);
    typedef double (*double_factorial_double)(const int &);
    m.def("__double_factorial", double_factorial_mp(&piranha::double_factorial), "Double factorial of non-negative integer argument.");
    m.def("__double_factorial", double_factorial_double(&piranha::double_factorial), "Double factorial of non-negative integer argument.");
    
    // Rising factorial.
    typedef double (*r_factorial_double)(const double &, const int &);
    typedef std::complex<double> (*r_factorial_complex_double)(const std::complex<double> &, const int &);
    typedef piranha::mp_integer (*r_factorial_mp_int)(const piranha::mp_integer &, const int &);
    typedef piranha::mp_rational (*r_factorial_mp_rat)(const piranha::mp_rational &, const int &);
    m.def("__r_factorial", r_factorial_complex_double(&piranha::r_factorial), "Rising factorial.");
    m.def("__r_factorial", r_factorial_double(&piranha::r_factorial), "Rising factorial.");
    m.def("__r_factorial", r_factorial_mp_int(&piranha::r_factorial), "Rising factorial.");
    m.def("__r_factorial", r_factorial_mp_rat(&piranha::r_factorial), "Rising factorial.");
    
    // Falling factorial.
    typedef double (*f_factorial_double)(const double &, const int &);
    typedef std::complex<double> (*f_factorial_complex_double)(const std::complex<double> &, const int &);
    typedef piranha::mp_integer (*f_factorial_mp_int)(const piranha::mp_integer &, const int &);
    typedef piranha::mp_rational (*f_factorial_mp_rat)(const piranha::mp_rational &, const int &);
    m.def("__f_factorial", f_factorial_complex_double(&piranha::f_factorial), "Falling factorial.");
    m.def("__f_factorial", f_factorial_double(&piranha::f_factorial), "Falling factorial.");
    m.def("__f_factorial", f_factorial_mp_int(&piranha::f_factorial), "Falling factorial.");
    m.def("__f_factorial", f_factorial_mp_rat(&piranha::f_factorial), "Falling factorial.");
    
    // Choose function.
    typedef double (*choose_double)(const int &, const int &);
    typedef double (*choose_double_double)(const double &, const int &);
    typedef std::complex<double> (*c_choose_double)(const std::complex<int> &, const int &);
    typedef std::complex<double> (*c_choose_double_double)(const std::complex<double> &, const int &);
    typedef piranha::mp_integer (*choose_z)(const piranha::mp_integer &, const int &);
    typedef piranha::mp_rational (*choose_q)(const piranha::mp_rational &, const int &);
    m.def("__choose", c_choose_double_double(&piranha::choose), "Binomial coefficient (complex double-precision).");
    m.def("__choose", c_choose_double(&piranha::choose), "Binomial coefficient (complex double-precision).");
    m.def("__choose", choose_double_double(&piranha::choose), "Binomial coefficient (double-precision).");
    m.def("__choose", choose_double(&piranha::choose), "Binomial coefficient (double-precision).");
    m.def("__choose", choose_z(&piranha::choose), "Binomial coefficient (multiprecision integer).");
    m.def("__choose", choose_q(&piranha::choose), "Binomial coefficient (multiprecision rational).");

    typedef double (*double_gamma)(double);
    m.def("gamma", double_gamma(&boost::math::tgamma<double>), "Gamma function.");

    typedef int (*cs_int)(const int &);
    typedef int (*cs_z)(const piranha::mp_integer &);
    m.def("cs_phase", cs_int(&piranha::cs_phase), "Condon-Shortley phase = (-1)**arg1.");
    m.def("cs_phase", cs_z(&piranha::cs_phase), "Condon-Shortley phase = (-1)**arg1.");
}
