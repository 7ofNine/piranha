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

#ifndef PIRANHA_TASS17_H
#define PIRANHA_TASS17_H

#include <complex>

#include "../bits/astro.h"
#include "../bits/prectest.h"

namespace piranha
{
#ifdef _PIRANHA_FORTRAN
  extern "C"
  {
    extern void posired_(const double &dj, const int &i, double xyz[3], double vxyz[3]);
  }
#else
  inline void posired_(const double &dj, const int &i, double xyz[3], double vxyz[3])
  {
    std::cout << "WARNING: FORTRAN support was not compiled in, returning 0." << std::endl;
    xyz[0]=xyz[1]=xyz[2]=vxyz[0]=vxyz[1]=vxyz[2]=0.;
  }
#endif

/// Class for TASS series.
  template <class Ps>
    class tass17
  {
      typedef std::complex<Ps> PsC;
    public:
/// Load into memory TASS series.
/**
 * The series are loaded in the same form seen in Alain Vienne's paper. They still are missing the
 * addition of the non-linear parts of the lambdas (i.e., the long period perturbations).
 * To correct the series use tass17::add_delta_lambdas.
 * If called more than once, the series will be reset to their default values.
 * @see tass17::add_delta_lambdas, to add long period perturbations to series.
 */
      static void load()
      {
        if (loaded)
        {
          std::cout << "Series already loaded, resetting them to the default values." << std::endl;
        }
        lambda4_=Ps("tass_lambda4.csv");
        lambda6_=Ps("tass_lambda6.csv");
        p6_=Ps("tass_p6.csv");
        z6_=PsC(Ps("tass_z6_real.csv"),Ps("tass_z6_imag.csv"));
//zeta6_=PsC(tass_zeta6); FIXME -> series is not available in electronic format yet.
        dlambda1_=Ps("tass_dlambda1.csv");
        dlambda2_=Ps("tass_dlambda2.csv");
        dlambda3_=Ps("tass_dlambda3.csv");
        dlambda4_=Ps("tass_dlambda4.csv");
        dlambda5_=Ps("tass_dlambda5.csv");
        dlambda6_=Ps("tass_dlambda6.csv");
        dlambda8_=Ps("tass_dlambda8.csv");
        loaded=true;
        has_deltas=false;
      }
/// Print to screen useful info about theory.
        static void status()
        {
          std::cout << "This is TASS version 1.7." << std::endl;
          std::cout << "-------------------------" << std::endl;
          std::cout << "Are series loaded?\t" << loaded << std::endl;
          std::cout << "Have deltas been added?\t" << has_deltas << std::endl;
        }
      static Ps r6()
      {
        std::pair<bool,psym_p> res=psymbol_manager::get_pointer("\\lambda_{o6}");
        if (!(res.first))
        {
          std::cout << "ERROR: no symbol named \\lambda_{o6} found, returning default series." << std::endl;
          return Ps();
        }
        const psym_p p=res.second;
        const double N=p->freq();
        Ps e6=e(z6());
        PsC complexp_M=eiM(lambda6(),z6(),e6);
        Ps cosE=astro::kep_cosE(e6,complexp_M,Ps::get_truncation());
        e6*=cosE;
        e6*=-1;
        e6+=1;
        e6*=a(N,p6());
        return e6;
      }
/// Calculate radius with Alain Vienne's FORTRAN routine.
/**
 * If FORTRAN support was not compiled in a warning message is displayed and the function returns 0.
 * @param[in] dj julian date.
 * @param[in] i satellite number.
 */
      static double vienne_r(const double &dj, int i)
      {
        double xyz[3], vxyz[3];
        posired_(dj,i,xyz,vxyz);
        return std::sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]+xyz[2]*xyz[2]);
      }
// Getters.
      static const double &m0()
      {
        return m0_;
      }
      static const double &m6()
      {
        return m6_;
      }
      static const Ps &lambda4()
      {
        return lambda4_;
      }
      static const Ps &lambda6()
      {
        return lambda6_;
      }
      static const Ps &p6()
      {
        return p6_;
      }
      static const PsC &z6()
      {
        return z6_;
      }
      static const PsC &zeta6()
      {
        return zeta6_;
      }
/// Add \f$ \delta\lambda_i \f$ to all series.
/**
 * Add the long period perturbations (\f$ \delta\lambda_i \f$ in Vienne's papers) to the series
 * of the theory.
 */
      static void add_delta_lambdas()
      {
        if (!loaded)
        {
          std::cout << "Please load series before adding deltas." << std::endl;
        }
        else if (!has_deltas)
        {
          add_delta_lambda(lambda4_);
          add_delta_lambda(lambda6_);
          add_delta_lambda(p6_);
          add_delta_lambda(z6_);
    //PsC tass17::zeta6_; FIXME!
          has_deltas=true;
        }
        else
        {
          std::cout << "Deltas have already been added, doing nothing." << std::endl;
        }
      }
/// Convert elliptic orbital element z into eccentricity e.
/**
 * This function simply calculates the absolute value of the complex iPsut series.
 * @param[in] z iPsut complex series for elliptic orbital element z.
 */
      static Ps e(const PsC &z)
      {
        return z.abs();
      }
/// Convert elliptic orbital element p into semi-major axis a.
/**
 * The returned value is expressed in the same time unit of input N to the power of \f$-\frac{2}{3}\f$.
 * Multiply by \f$ \sqrt[3]{\mu} \f$ to get a metric quantity.
 * \f$ \mu \f$ is the sum of the
 * gravitational parameters of Saturn and of the satellite whose semi-major axis we are calculating.
 * @param[in] N double precision mean motion.
 * @param[in] p input real series for the elliptic orbital element p.
 */
      static Ps a(const double &N, const Ps &p)
      {
        Ps a=p;
        a+=1;
        a*=N;
        a=a.pow(-2./3.);
        return a;
      }
/// Find complex exponential of mean mean motion M.
      static PsC eiM(const Ps &lambda, const PsC &z, const Ps &e)
      {
        PsC retval=lambda.complexp();
        retval*=z.conj();
        retval*=e.pow(-1.);
        return retval;
      }
    private:
      template <class T>
        static void add_delta_lambda(T &p)
      {
        p.add_ps_to_arg("\\lambda_{o1}",dlambda1_);
        p.add_ps_to_arg("\\lambda_{o2}",dlambda2_);
        p.add_ps_to_arg("\\lambda_{o3}",dlambda3_);
        p.add_ps_to_arg("\\lambda_{o4}",dlambda4_);
        p.add_ps_to_arg("\\lambda_{o5}",dlambda5_);
        p.add_ps_to_arg("\\lambda_{o6}",dlambda6_);
        p.add_ps_to_arg("\\lambda_{o8}",dlambda8_);
      }
// Data members.
    private:
// Tass physical parameters.
      static const double m0_;
      static const double m6_;
      static Ps lambda4_;
      static Ps lambda6_;
      static Ps p6_;
      static PsC z6_;
      static PsC zeta6_;
      static Ps dlambda1_;
      static Ps dlambda2_;
      static Ps dlambda3_;
      static Ps dlambda4_;
      static Ps dlambda5_;
      static Ps dlambda6_;
      static Ps dlambda8_;
// Flag to see if series were already loaded.
      static bool loaded;
// Flag to see if deltas have been added
      static bool has_deltas;
  };

  template <class Ps>
    class tc_vienne_r6:public base_tc<Ps>
  {
    public:
// b_type stands for "benchmarked type"
      typedef Ps b_type;
      typedef typename base_tc<Ps>::eval_type eval_type;
      tc_vienne_r6(const Ps &b, const double &t1, const double &t2, const size_t &ntot)
      {
        base_tc<Ps>::init(t1,t2,ntot,b);
      }
    private:
// t is expressed in Julian years from J1980.0.
      virtual eval_type eval_hs_computed(const double &t) const
      {
        return tass17<Ps>::vienne_r((astro::J1980dot0()+t*astro::JD_per_JY()),6)*
          astro::AU();
      }
  };

// Initialization of tass17's static variables
  template <class Ps>
    bool tass17<Ps>::loaded=false;
  template <class Ps>
    bool tass17<Ps>::has_deltas=false;
  template <class Ps>
    const double tass17<Ps>::m0_=3498.790;
  template <class Ps>
    const double tass17<Ps>::m6_=0.4225863977890E+04;
  template <class Ps>
    Ps tass17<Ps>::lambda4_;
  template <class Ps>
    Ps tass17<Ps>::lambda6_;
  template <class Ps>
    Ps tass17<Ps>::p6_;
  template <class Ps>
    typename tass17<Ps>::PsC tass17<Ps>::z6_;
  template <class Ps>
    typename tass17<Ps>::PsC tass17<Ps>::zeta6_;
  template <class Ps>
    Ps tass17<Ps>::dlambda1_;
  template <class Ps>
    Ps tass17<Ps>::dlambda2_;
  template <class Ps>
    Ps tass17<Ps>::dlambda3_;
  template <class Ps>
    Ps tass17<Ps>::dlambda4_;
  template <class Ps>
    Ps tass17<Ps>::dlambda5_;
  template <class Ps>
    Ps tass17<Ps>::dlambda6_;
  template <class Ps>
    Ps tass17<Ps>::dlambda8_;
}
#endif
