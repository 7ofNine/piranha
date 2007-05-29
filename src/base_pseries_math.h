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

#ifndef PIRANHA_BASE_PSERIES_MATH_H
#define PIRANHA_BASE_PSERIES_MATH_H

namespace piranha
  {
  // Assignment operator
  // -------------------
  template <class Cf,class Trig,template <class,class> class I>
  inline void base_pseries<Cf,Trig,I>::basic_assignment(const
      base_pseries<Cf,Trig,I> &ps2)
  {
    if (this==&ps2)
      {
        return;
      }
    set_=ps2.set_;
    norm_=ps2.norm_;
    cf_s_vec_=ps2.cf_s_vec_;
    trig_s_vec_=ps2.trig_s_vec_;
    lin_args_=ps2.lin_args_;
    std::cout << "Assignment operator!" << std::endl;
  }

  /************************/
  /* Low level operations */
  /************************/

  // Base merge operator
  // -------------------
  template <class Cf,class Trig,template <class,class> class I>
  template <class Cf2,class Trig2,template <class,class> class I2>
  inline void base_pseries<Cf,Trig,I>::alg_sum_lin_args(const base_pseries<Cf2,Trig2,I2> &ps2,
      bool sign)
  {
    vector_mult_t tmp(trig_width());
    if (sign)
      {
        vec_add(lin_args_,ps2.lin_args(),ps2.lin_args(),lin_args_,tmp);
      }
    else
      {
        vec_sub(lin_args_,ps2.lin_args(),ps2.lin_args(),lin_args_,tmp);
      }
    lin_args_=tmp;
  }

  template <class Cf,class Trig,template <class,class> class I>
  template <class Cf2,class Trig2,template <class,class> class I2>
  inline void base_pseries<Cf,Trig,I>::merge_with(const base_pseries<Cf2,Trig2,I2> &ps2, bool sign)
  {
    if ((void *)&ps2==(void *)this)
      {
        if (sign)
          {
            base_pseries tmp_ps(*this);
            tmp_ps.merge_with(ps2,sign);
            swap(tmp_ps);
          }
        else
          {
            base_pseries tmp_ps;
            tmp_ps.merge_args(*this);
            tmp_ps.lin_args_=lin_args();
            swap(tmp_ps);
          }
        return;
      }
    // Check that trig_args are compatible
    if (!merge_args(ps2))
      {
        std::cout << "trig_args are not compatible, returning self." << std::endl;
        std::exit(1);
        return;
      }
    // Sum/sub lin_args
    alg_sum_lin_args(ps2,sign);
    // Use hint, since as we add terms we have an idea of where they are going to be placed
    it_s_index it_hint=s_index().end();
    arg_manager::arg_assigner aa(&cf_s_vec_,&trig_s_vec_);
    // NOTE: At this point this' size is greater or equal to ps2'
    for (typename base_pseries<Cf2,Trig2,I2>::it_h_index it=ps2.h_index().begin();
         it!=ps2.h_index().end();++it)
      {
        it_hint=insert(*it,sign,&it_hint);
      }
  }


  // Merge with a generic entity - NOT with another series
  // -----------------------------------------------------
  template <class Cf,class Trig,template <class,class> class I>
  template <class T>
  inline void base_pseries<Cf,Trig,I>::generic_merge(const T &x)
  {
    // Build a series from x
    base_pseries tmp=base_pseries(cf_type(x));
    // Merge with this
    merge_with(tmp);
  }


  // Low-level mutiplication of terms
  // --------------------------------
  template <class Cf,class Trig,template <class,class> class I>
  template <class Cf2,class Trig2,template <class,class> class I2>
  inline void base_pseries<Cf,Trig,I>::mult_terms(const base_pseries<Cf2,Trig2,I2> &ps2,
      base_pseries &retval, const double &Delta)
  {
    const double Delta_threshold=Delta/(2*length()*ps2.length());
    unsigned long int n=0;
    // NOTE: at this point retval's width() is greater or equal to _both_ this
    // and ps2. It's the max width indeed.
    p_assert(math::max(trig_width(),ps2.trig_width())==retval.trig_width());
    term_type tmp1, tmp2;
    boost::tuple<term_type &, term_type &> term_pair(tmp1,tmp2);
    const it_s_index it1_f=s_index().end();
    const typename base_pseries<Cf2,Trig2,I2>::it_s_index it2_f=ps2.s_index().end();
    typename base_pseries<Cf2,Trig2,I2>::it_s_index it2;
    it_s_index it1, it_hint=retval.s_index().end();
    for (it1=s_index().begin();it1!=it1_f;++it1)
      {
        it2=ps2.s_index().begin();
        if ((it1->norm(cf_s_vec_)*it2->norm(cf_s_vec_))/2<Delta_threshold)
          {
            break;
          }
        for (;it2!=it2_f;++it2)
          {
            // We are going to calculate a term's norm twice... We need to profile
            // this at a later stage and see if it is worth to store the norm inside
            // the term.
            if ((it1->norm(cf_s_vec_)*it2->norm(cf_s_vec_))/2<Delta_threshold)
              {
                break;
              }
            it1->mult_by(*it2,term_pair);
            // Before insertion we change the sign of trigonometric parts if necessary.
            // This way we won't do a copy inside insertion function.
            if (term_pair.template get
                  <0>().trig_args().sign()<0)
                {
                  term_pair.template get<0>().invert_trig_args();
                }
            if (term_pair.template get
                  <1>().trig_args().sign()<0)
                {
                  term_pair.template get<1>().invert_trig_args();
                }
            it_hint=retval.insert(term_pair.template get
                                    <0>(),true,&it_hint);
            it_hint=retval.insert(term_pair.template get
                                    <1>(),true,&it_hint);
            ++n;
          }
      }
    //retval.cumulative_crop(Delta);
    std::cout << "w/o trunc=" << length()*ps2.length() << "\tw/ trunc=" << n << std::endl;
    std::cout << "Out length=" << retval.length() << std::endl;
  }



  // Basic multiplication
  // --------------------
  template <class Cf,class Trig,template <class,class> class I>
  template <class T>
  inline void base_pseries<Cf,Trig,I>::basic_ps_mult(const
      T &ps2)
  {
    base_pseries retval;
    // OPTIMIZE: optimize in the case one series is a c value
    // If one length is zero do not do anything
    if (length()!=0 && ps2.length()!=0)
      {
        if (!is_zero_vec(lin_args_)||!is_zero_vec(ps2.lin_args()))
          {
            std::cout << "Non-zero linargs!" << std::endl;
            std::exit(1);
          }
        p_assert(retval.merge_args(*this));
        if (!retval.merge_args(ps2))
          {
            std::cout << "mult trig_argss are not compatible, returning self." << std::endl;
            std::exit(1);
            return;
          }
        const double Delta=norm()*ps2.norm()*settings_manager::prec();
        arg_manager::arg_assigner aa(&cf_s_vec_,&trig_s_vec_);
        mult_terms(ps2,retval,Delta);
      }
    else
      {
        std::cout << "Zero stuff" << std::endl;
      }
    swap(retval);
  }


  // Multiplication by a generic entity
  // ----------------------------------
  template <class Cf,class Trig,template <class,class> class I>
  template <class T>
  inline void base_pseries<Cf,Trig,I>::generic_mult(const T &c)
  {
    if (length()==0)
      {
        return;
      }
    if (!is_zero_vec(lin_args_))
      {
        std::cout << "Non-zero linargs in *= T!" << std::endl;
        std::exit(1);
      }
    base_pseries tmp_ps;
    tmp_ps.merge_args(*this);
    term_type tmp_term;
    it_s_index it_hint=tmp_ps.s_index().end();
    arg_manager::arg_assigner aa(&cf_s_vec_,&trig_s_vec_);
    for (it_s_index it=s_index().begin();it!=s_index().end();++it)
      {
        tmp_term=*it;
        tmp_term.c()*=c;
        it_hint=tmp_ps.insert(tmp_term,true,&it_hint);
      }
    swap(tmp_ps);
  }


  // Multiplication by an integer
  // ----------------------------
  template <class Cf,class Trig,template <class,class> class I>
  inline void base_pseries<Cf,Trig,I>::mult_by_int(int n)
  {
    const vector_mult_t old_lin_args=lin_args_;
    unsigned int j;
    // Zero the linargs, otherwise the generic *= operator complains
    for (j=0;j<lin_args_.size();++j)
      {
        lin_args_[j]=0;
      }
    // Now perform the generic multiplication
    generic_mult(n);
    // Multiply the old linargs and restore them
    const unsigned int w=lin_args_.size();
    for (j=0;j<w;++j)
      {
        lin_args_[j]=old_lin_args[j]*n;
      }
  }


  /************************/
  /* Advanced operations  */
  /************************/

  // Real power.

  // Calculate limit of the development given the desired error and the power
  // Requirements: error >= 0.

#define __pow_hard_limit 20
  template <class Cf,class Trig,template <class,class> class I>
  inline unsigned int base_pseries<Cf,Trig,I>::pow_limit(const double &error,
      const double &power) const
    {
      p_assert(error>=0);
      unsigned int retval=0;
      const double a=s_index().begin()->norm(cf_s_vec_), absx=norm()-a,
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
  template <class Cf,class Trig,template <class,class> class I>
  inline void base_pseries<Cf,Trig,I>::basic_pow(const double &power)
  {
    if (length()==0)
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
        return;
      }
    if (s_index().begin()->flavour()==false ||
        !(s_index().begin()->trig_args().is_zero()))
      {
        std::cout << "ERROR! series' top term is not suitable for real power." << std::endl;
        std::exit(1);
      }
    // NOTICE: what does it mean to evaluate here for symbolic coefficients? To be really effective
    // symbolic coefficients should not have any time dependency. Otherwise this is just an approximation.
    // Need to think about this, but it is not essential until symbolic coefficients are introduced.
    const cf_type &a=s_index().begin()->c();
    if (a.t_eval(0.,cf_s_vec_)<0)
      {
        std::cout << "ERROR! I want a positive evaluation for the greatest term's coefficient." << std::endl;
        std::exit(1);
      }
    // Top term must be greater than half of the series' norm.
    if (2*a.norm(cf_s_vec_)<=norm())
      {
        std::cout << "ERROR! series' top term is not big enough for negative power." << std::endl;
        std::exit(1);
      }
    // NOTICE: Hard coded binomial expansion error to 1/10 of desired precision.
    const double error=.1*std::pow(norm(),power)*settings_manager::prec();
    const unsigned int limit_index=pow_limit(error,power);
    base_pseries retval, x(*this), tmp(cf_type(1.));
    x.term_erase(x.s_index().begin());
    for (unsigned int i=0;i<=limit_index;++i)
      {
        {
          base_pseries tmp2(tmp);
          tmp2*=math::choose(power,i);
          tmp2*=a.pow(power-i);
          retval+=tmp2;
        }
        tmp*=x;
      }
    swap(retval);
  }
}

#endif
