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

#include "../astro.h"
#include "../lnp.h"
#include "../lnpc.h"
#include "../prectest.h"

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
  class tass17
    {
    public:
      static void load();
      static void status();
      static lnp e(const lnpc &);
      static lnp a(const double &, const lnp &);
      static lnpc eiM(const lnp &, const lnpc &, const lnp &);
      static lnp r6();
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
      static const lnp &lambda4()
      {
        return lambda4_;
      }
      static const lnp &lambda6()
      {
        return lambda6_;
      }
      static const lnp &p6()
      {
        return p6_;
      }
      static const lnpc &z6()
      {
        return z6_;
      }
      static const lnpc &zeta6()
      {
        return zeta6_;
      }
      // Add delta lambdas to all series.
      static void add_delta_lambdas();
    private:
      template <class T>
      static void add_delta_lambda(T &);
      // Data members.
    private:
      // Tass physical parameters.
      static const double m0_;
      static const double m6_;

      static lnp lambda4_;
      static lnp lambda6_;

      static lnp p6_;

      static lnpc z6_;

      static lnpc zeta6_;

      static lnp dlambda1_;
      static lnp dlambda2_;
      static lnp dlambda3_;
      static lnp dlambda4_;
      static lnp dlambda5_;
      static lnp dlambda6_;
      static lnp dlambda8_;
      // Flag to see if series were already loaded.
      static bool loaded;
      // Flag to see if deltas have been added
      static bool has_deltas;
    };

  class tc_vienne_r6:public base_tc<lnp>
    {
    public:
      // b_type stands for "benchmarked type"
      typedef lnp b_type;
      typedef base_tc<lnp>::eval_type eval_type;
      tc_vienne_r6(const lnp &b, const double &t1, const double &t2, const size_t &ntot)
      {
        base_tc<lnp>::init(t1,t2,ntot,b);
      }
    private:
      // t is expressed in Julian years from J1980.0.
      virtual eval_type eval_hs_computed(const double &t) const
        {
          return tass17::vienne_r((astro::J1980dot0()+t*astro::JD_per_JY()),6)*
                 astro::AU();
        }
    };


  inline lnp tass17::r6()
  {
    psym_p p=psymbol_manager::get_pointer("\\lambda_{o6}");
    if (p==psymbol_manager::end())
      {
        std::cout << "ERROR: no symbol named \\lambda_{o6} found, returning defaul series." << std::endl;
        return lnp();
      }
    const double N=p->freq();
    lnp e6=e(z6());
    lnpc complexp_M=eiM(lambda6(),z6(),e6);
    lnp cosE=astro::kep_cosE(e6,complexp_M,settings_manager::prec());
    e6*=cosE;
    e6*=-1;
    e6+=1;
    e6*=a(N,p6());
    return e6;
  }


  /// Load into memory TASS series.
  /**
   * The series are loaded in the same form seen in Alain Vienne's paper. They still are missing the
   * addition of the non-linear parts of the lambdas (i.e., the long period perturbations).
   * To correct the series use tass17::add_delta_lambdas.
   * If called more than once, the series will be reset to their default values.
   * @see tass17::add_delta_lambdas, to add long period perturbations to series.
   */
  inline void tass17::load()
  {
    if (loaded)
      {
        std::cout << "Series already loaded, resetting them to the default values." << std::endl;
      }

    lambda4_=lnp("tass_lambda4.csv");
    lambda6_=lnp("tass_lambda6.csv");

    p6_=lnp("tass_p6.csv");

    z6_=lnpc("tass_z6_real.csv","tass_z6_imag.csv");

    //zeta6_=lnpc(tass_zeta6); FIXME -> series is not available in electronic format yet.

    dlambda1_=lnp("tass_dlambda1.csv");
    dlambda2_=lnp("tass_dlambda2.csv");
    dlambda3_=lnp("tass_dlambda3.csv");
    dlambda4_=lnp("tass_dlambda4.csv");
    dlambda5_=lnp("tass_dlambda5.csv");
    dlambda6_=lnp("tass_dlambda6.csv");
    dlambda8_=lnp("tass_dlambda8.csv");

    loaded=true;
    has_deltas=false;
  }


  /// Print to screen useful info about theory.
  inline void tass17::status()
  {
    std::cout << "This is TASS version 1.7." << std::endl;
    std::cout << "-------------------------" << std::endl;
    std::cout << "Are series loaded?\t" << loaded << std::endl;
    std::cout << "Have deltas been added?\t" << has_deltas << std::endl;
  }


  template <class T>
  inline void tass17::add_delta_lambda(T &p)
  {
    p.add_ps_to_arg("\\lambda_{o1}",dlambda1_);
    p.add_ps_to_arg("\\lambda_{o2}",dlambda2_);
    p.add_ps_to_arg("\\lambda_{o3}",dlambda3_);
    p.add_ps_to_arg("\\lambda_{o4}",dlambda4_);
    p.add_ps_to_arg("\\lambda_{o5}",dlambda5_);
    p.add_ps_to_arg("\\lambda_{o6}",dlambda6_);
    p.add_ps_to_arg("\\lambda_{o8}",dlambda8_);
  }


  /// Add \f$ \delta\lambda_i \f$ to all series.
  /**
   * Add the long period perturbations (\f$ \delta\lambda_i \f$ in Vienne's papers) to the series
   * of the theory.
   */
  inline void tass17::add_delta_lambdas()
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

        //lnpc tass17::zeta6_; FIXME!

        has_deltas=true;
      }
    else
      {
        std::cout << "Deltas have already been added, doing nothing." << std::endl;
      }
  }


  /// Convert elliptic orbital element z into eccentricity e.
  /**
   * This function simply calculates the absolute value of the complex ilnput series.
   * @param[in] z ilnput complex series for elliptic orbital element z.
   */
  inline lnp tass17::e(const lnpc &z)
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
  inline lnp tass17::a(const double &N, const lnp &p)
  {
    lnp a=p;
    a+=1.;
    a*=N;
    a=a.pow(-2./3.);
    return a;
  }


  /// Find complex exponential of mean mean motion M.
  inline lnpc tass17::eiM(const lnp &lambda, const lnpc &z, const lnp &e)
  {
    lnpc retval=lambda.complexp();
    retval*=z.conj();
    retval*=e.pow(-1.);
    return retval;
  }
}

#endif
