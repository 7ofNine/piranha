/***************************************************************************
 *   Copyright (C) 2007 by Francesco Biscani   *
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

#ifndef PIRANHA_ASTRO_H
#define PIRANHA_ASTRO_H

#include "math.h"
#include "settings_manager.h"

namespace piranha
{
/// Useful astronomical functions and constants.
  class astro
  {
    public:
/// Transform given Julian Date into time unit suitable for use in the ELP2000 theory.
/**
 * ELP2000 time is measured in Julian centuries from J2000.0.
 * @param[in] jd double: Julian Date to be converted.
 */
      static double JD_to_elp2000(const double &jd)
      {
        return ((jd-J2000dot0_)/36525.);
      }
/// Kepler equation solver for \f$ \cos E \f$.
/**
 * Find the cosine of eccentric anomaly by successive approximations,
 * given eccentricity and complex exponential of mean anomaly.
 * @param[in] e: eccentricity.
 * @param[in] eiM: complex exponential of mean anomaly.
 * @param[in] prec: double precision desired for the solution, relative to \f$ E \f$.
 */
      template <class T>
        static T kep_cosE(const T &e,
        const std::complex<T> &eiM, const double &prec)
      {
        if (prec<0 || prec>1 || std::abs(prec)<settings_manager::numerical_zero())
        {
          std::cout << "Error: invalid precision requested in kep_cosE" << std::endl;
          std::exit(1);
        }
        if (math::norm(e)-settings_manager::numerical_zero()>=1)
        {
          std::cout << "Error: invalid eccentricity requested in kep_cosE" << std::endl;
          std::exit(1);
        }
        const T cosM=eiM.real(), sinM=eiM.imag();
        int n=(int)std::ceil((std::log10(prec)/std::log10(math::norm(e)))-1);
        p_assert(n>=0);
        T phi(int(0));
        std::complex<T> eiphi=math::complexp(phi);
        for (int i=0;i<=n;++i)
        {
          phi=sinM;
          phi*=eiphi.real();
          phi+=cosM*eiphi.imag();
          phi*=e;
          eiphi=math::complexp(phi);
        }
        return (cosM*eiphi.real()-sinM*eiphi.imag());
      }
/// Transform spherical coordinates into x-coordinate.
      static double sph_to_x(const double &r, const double &col, const double &lon)
      {
        return (r*std::cos(lon)*std::sin(col));
      }
/// Transform spherical coordinates into y-coordinate.
      static double sph_to_y(const double &r, const double &col, const double &lon)
      {
        return (r*std::sin(lon)*std::sin(col));
      }
/// Transform spherical coordinates into z-coordinate.
      static double sph_to_z(const double &r, const double &col, const double &)
      {
        return (r*std::cos(col));
      }
// Getters for constants.
/// Get universal gravitational constant.
      static const double &G()
      {
        return G_;
      }
/// Get gaussian gravitational constant.
      static const double &k()
      {
        return k_;
      }
/// Get Julian Date for J2000.0.
      static const double &J2000dot0()
      {
        return J2000dot0_;
      }
/// Get Julian Date for J1980.0.
      static const double &J1980dot0()
      {
        return J1980dot0_;
      }
/// Get number of Julian Days per Julian Year.
      static const double &JD_per_JY()
      {
        return JD_per_JY_;
      }
/// Get number of seconds per Julian Year.
      static const double &seconds_per_JY()
      {
        return seconds_per_JY_;
      }
// FIXME: is this really at J2000.0?
/// Get obliquity of the ecliptic at J2000.0.
      static const double &eps_0()
      {
        return eps_0_;
      }
/// Get the Astronomical Unit in meters.
      static const double &AU()
      {
        return AU_;
      }
    private:
      static const double G_;
      static const double k_;
      static const double J2000dot0_;
      static const double J1980dot0_;
      static const double JD_per_JY_;
      static const double seconds_per_JY_;
      static const double eps_0_;
      static const double AU_;
  };

/*inline boost::python::numeric::array sph_to_cart(const double &r, const double &col, const double &lon)
{
    return boost::python::numeric::array(1.);
}*/
}
#endif
