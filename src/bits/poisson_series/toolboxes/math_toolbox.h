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

#ifndef PIRANHA_MATH_TOOLBOX_H
#define PIRANHA_MATH_TOOLBOX_H

namespace piranha
{
  template <class Derived>
    class math_toolbox
  {
    public:
      /// Real power.
      /**
       * Calculate the power of a series. The power can be any real number. It employs the generalized
       * Newton binomial to expand the real power into a series of natural powers. The series' maximum
       * term must have certain features: it must be a cosine with all argument indices set to zero, and
       * the evaluation of its coefficient must be positive and greater than the sum of the norms of all
       * remaining terms. Basically it is as effective as a Taylor expansion.
       *
       * It may be possible in the future to extend to negative coefficients, but in that case the output
       * will have to be a complex series.
       * @param[in] power real power the series will be raised to.
       */
      Derived pow(const double &power) const
      {
        typedef typename Derived::ancestor::cf_type cf_type;
        const Derived *derived_cast=static_cast<Derived const *>(this);
        if (derived_cast->length()==0)
        {
          if (std::abs(power)<settings_manager::numerical_zero())
          {
            std::cout << "WARNING: won't raise nil power to zero, returning self." << std::endl;
            std::exit(1);
          }
          else if (power<0)
          {
            std::cout << "ERROR! won't neg power a nil series." << std::endl;
            std::exit(1);
          }
          return Derived();
        }
        if (!(derived_cast->g_s_index().begin()->trig().flavour()) or
          !(derived_cast->g_s_index().begin()->trig().is_zero()))
        {
          std::cout << "ERROR! series' top term is not suitable for real power." << std::endl;
          std::exit(1);
        }
        // NOTICE: what does it mean to evaluate here for symbolic coefficients? To be really effective
        // symbolic coefficients should not have any time dependency. Otherwise this is just an approximation.
        // Need to think about this, but it is not essential until symbolic coefficients are introduced.
        const cf_type &a=derived_cast->g_s_index().begin()->cf();
        if (a.t_eval(0.,derived_cast->arguments().template get<0>())<0)
        {
          std::cout << "ERROR! I want a positive evaluation for the greatest term's coefficient." << std::endl;
          std::exit(1);
        }
        // Top term must be greater than half of the series' norm.
        if (2*a.norm(derived_cast->arguments().template get<0>()) <= derived_cast->g_norm())
        {
          std::cout << "ERROR! series' top term is not big enough for negative power." << std::endl;
          std::exit(1);
        }
        // NOTICE: Hard coded binomial expansion error to 1/10 of desired precision.
        const double error=.1*std::pow(derived_cast->g_norm(),power)*
          derived_cast->get_truncation();
        const unsigned int limit_index=pow_limit(error,power);
        Derived retval, x(*derived_cast), tmp(1);
        x.term_erase(x.g_s_index().begin());
        for (unsigned int i=0;i<=limit_index;++i)
        {
          {
            Derived tmp2(tmp);
            tmp2.mult_by(math::choose(power,i));
            tmp2.mult_by(Derived(a.pow(power-i),*derived_cast));
            retval+=tmp2;
          }
          tmp*=x;
        }
        return retval;
      }
    private:
      /// Real power limit.
      /**
       * Calculate limit of the development given the desired error and the power. Requirements: error >= 0.
       */
#define __pow_hard_limit 20
      unsigned int pow_limit(const double &error, const double &power) const
      {
        p_assert(error>=0);
        unsigned int retval=0;
        const Derived *derived_cast=static_cast<Derived const *>(this);
        const double a=derived_cast->g_s_index().begin()->
          cf().norm(derived_cast->arguments().template get<0>()),
          absx=derived_cast->g_norm()-a,
          exactM=std::pow(a+absx,power), exactm=std::pow(a-absx,power), absratio=absx/a;
        p_assert(a>0);
        p_assert(a>absx);
        double DeltaM=exactM, Deltam=exactm;
        do
        {
          Deltam-=math::choose(power,retval)*math::natural_pow(retval,-absratio)*std::pow(a,power);
          DeltaM-=math::choose(power,retval)*math::natural_pow(retval,absratio)*std::pow(a,power);
          std::cout << "Deltam=" << Deltam << '\n';
          std::cout << "DeltaM=" << DeltaM << '\n';
          if (std::max(std::abs(DeltaM),std::abs(Deltam))<error)
          {
            break;
          }
          ++retval;
        }
        while (retval<__pow_hard_limit);
        return retval;
      }
#undef __pow_hard_limit
  };
}
#endif
