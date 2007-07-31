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

#ifndef PIRANHA_MATH_H
#define PIRANHA_MATH_H

#include <complex>
#include <iostream>
#include <vector>

#include "p_assert.h"
#include "platform_switches.h"                    // For NAN detection.
#include "settings_manager.h"

// Sign function - FIXME: turn into a template?
#define __sgn(x) (x==0.0 ? 0.0 : (x > 0.0 ? 1.0 : -1.0))

namespace piranha
{
/// Class that groups together useful mathematical functions.
  class math
  {
    public:
// Operations on vectors
// Add/subtract two differently sized containers, whose members can be accessed through the
// "[]" operator.
// The result goes into retval.
// Pre-requisites: v1.size() >= v2.size(), retval.size()==v1.size().
      template <class T>
        static void vec_add(const T &v1, const T &v2, T &retval)
      {
        const size_t w1=v1.size(), w2=v2.size();
        p_assert(w1>=w2);
        p_assert(w1==retval.size());
        size_t j;
        for (j=0;j<w2;++j)
        {
          retval[j]=v1[j]+v2[j];
        }
        for (;j<w1;++j)
        {
          retval[j]=v1[j];
        }
      }
// Pre-requisites: v1.size() >= v2.size(), retval.size()==v1.size().
      template <class T>
        static void vec_sub(const T &v1, const T &v2, T &retval)
      {
        const size_t w1=v1.size(), w2=v2.size();
        p_assert(w1>=w2);
        p_assert(w1==retval.size());
        size_t j;
        for (j=0;j<w2;++j)
        {
          retval[j]=v1[j]-v2[j];
        }
        for (;j<w1;++j)
        {
          retval[j]=v1[j];
        }
      }
/// Check whether a vector is filled with zero elements.
      template <class T>
        static bool is_zero_vec(const T &v)
      {
        const size_t w=v.size();
        for (size_t j=0;j<w;++j)
        {
          if (v[j]!=0)
          {
            return false;
          }
        }
        return true;
      }
/// Sign function.
      static short int sgn(int x)
      {
        if (x<0)
        {
          return -1;
        }
        else
        {
          return 1;
        }
      }
/// Kroenecker Delta function.
      template <class T, class U>
        static int Kdelta(const T &n1, const U &n2)
      {
        return (n1==n2);
      }
/// Find maximum value.
      template <class T>
        static const T &max(const T &a, const T &b)
      {
        if (a>b)
        {
          return a;
        }
        else
        {
          return b;
        }
      }
/// Find minimum value.
      template <class T>
        static const T &min(const T &a, const T &b)
      {
        if (a<b)
        {
          return a;
        }
        else
        {
          return b;
        }
      }
/// Check whether an integer is odd.
      static bool is_odd(const long int &n)
      {
        return (n&1);
      }
/// Return norm, generic version.
      template <class T>
        static double norm(const T &x)
      {
        return x.g_norm();
      }
/// Condon-Shortley phase.
      static short int cs_phase(const long int &n)
      {
        if (is_odd(n))
        {
          return -1;
        }
        else
        {
          return 1;
        }
      }
/// Generalized binomial coefficient.
      template <class T>
        static T choose(const T &z, const long int &k)
      {
// k cannot be negative.
        if (k<0)
        {
          return T(0);
        }
        T retval(1);
        for (long int i=1;i<=k;++i)
        {
          retval*=(z-T(k)+T(i));
          retval/=T(i);
        }
        return retval;
      }
/// Factorial.
      static double factorial(int n)
      {
        return ll_factorial<double>(n);
      }
/// Generic factorial.
      template <class T>
        static T generic_factorial(int n)
      {
        return ll_factorial<T>(n);
      }
/// Double factorial.
      static double dbl_factorial(int ext_n)
      {
        double k=1;
        unsigned int i=1;
        if (ext_n>=0)
        {
          unsigned int n=ext_n;
          for (i=1;i<n;i=(i+2))
          {
            k*=(i+2);
          }
        }
        else if (ext_n==-1)
          k=1;
        else
        {
          unsigned int n=-ext_n-2;
          for (i=1;i<n;i=(i+2))
          {
            k*=(i+2);
          }
          k=-1./k;
        }
        return k;
      }
/// Bessel J series limit.
/**
 * Identifies the number of iterations needed to compute a Bessel function of the first kind
 * through the power series definition in order to reach a precision of the same order of magnitude of
 * Piranha's general relative precision.
 * @param[in] n_, integer order of the Bessel function.
 * @param[in] x, double argument of the Bessel function.
 */
#define __max_bessel_iter (20)
      static unsigned int besselJ_series_limit(int n_, const double &x)
      {
        if (std::abs(x)<settings_manager::numerical_zero())
        {
          return 0;
        }
        unsigned int n;
        short int sign_mult;
        double fact;
        if (n_<0)
        {
          n=-n_;
          sign_mult=cs_phase(n);
        }
        else
        {
          n=n_;
          sign_mult=1;
        }
        if (n==0)
        {
          fact=1;
        }
        else
        {
          fact=generic_factorial<double>(n);
        }
        double half_x=(x/2);
        double tmp=1, a_jm1=1, a_j, half_x_pow2=natural_pow(2,half_x);
        unsigned int j=0;
        double tmp_besselJ=besselJ(n_,x), tmp_pow=natural_pow(n,half_x);
        double target=0;
        if (std::abs(tmp_besselJ)<settings_manager::numerical_zero())
        {
//std::cout << "Numerical limits of double reached, returning 0." << std::endl;
        }
        else if (std::abs(tmp_pow)<settings_manager::numerical_zero())
        {
// We have problems here since we will be dividing by tmp_pow. Exit with j=0.
// TODO: add precision warning.
//std::cout << "Numerical limits of natural_pow reached, returning 0." << std::endl;
        }
        else
        {
// Target is multiplied that way because the retval is modified after the "for" cycle,
// so we must "anticipate" the modification to test whether we are ok with the precision.
          target=(tmp_besselJ*fact)/tmp_pow*sign_mult;
          const double target_precision=target*settings_manager::prec();
          do
          {
//std::cout << "tmp is: " << tmp << '\n';
            ++j;
            (a_jm1/=(j*(n+j)))*=-1;
            a_j=a_jm1;
            a_j*=half_x_pow2;
//std::cout << "a_j is: " << a_j << '\n';
            if (std::abs(a_j)<settings_manager::numerical_zero())
            {
// We cannot improve precision any more, since a_j went to zero. Exit the cycle.
// TODO: add precision warning.
//std::cout << "a_j went zero, returning current j\n";
              break;
            }
            tmp+=a_j;
            a_jm1=a_j;
          }
          while (std::abs(tmp-target)>target_precision &&
            j <= __max_bessel_iter);
        }
        if (j>__max_bessel_iter)
        {
          std::cout << "OH NOES!\n";
          std::exit(1);
        }
        return j;
      }
#undef __max_bessel_iter
/// Bessel function of the first kind. Power series implementation.
/**
 * Uses the series definition for Bessel functions. To be used with symbolic types,
 * for numerical types it is much better to use math::besselJ.
 * @param[in] n_, integer order of the Bessel function.
 * @param[in] x, double argument of the Bessel function.
 * @param[in] iterations, unsigned int number of iterations.
 */
      template <class T, class Integer>
        static T pow_besselJ(int n_, const T &x, unsigned int iterations)
      {
        unsigned int n;
        int sign_mult;
        Integer fact;
        if (n_<0)
        {
          n=-n_;
          sign_mult=cs_phase(-n_);
        }
        else
        {
          n=n_;
          sign_mult=1;
        }
        if (n==0)
        {
          fact=1;
        }
        else
        {
          fact=ll_factorial<Integer>(n);
        }
        T half_x=x;
        half_x/=2;
        T retval=T(1), a_jm1=T(1), a_j, half_x_pow2=T(natural_pow(2,half_x));
        for (unsigned int j=1;j<=iterations;++j)
        {
          (a_jm1/=(int)(j*(n+j)))*=-1;
          a_j=a_jm1;
          a_j*=half_x_pow2;
          retval+=a_j;
          a_jm1=a_j;
        }
        retval*=natural_pow(n,half_x);
        retval/=fact;
        retval*=sign_mult;
        return retval;
      }
/// Bessel function of the first kind (C standard library implementation).
      static double besselJ(int n, const double &x)
      {
        double retval=jn(n,x);
// Check needed apparently under MinGW/GCC.
        if (__ISNAN(retval))
        {
          std::cout << "I don't believe it!\n";
          retval=0.;
        }
        return retval;
      }
/// Legendre function of the first kind - Pnm(cos(theta)).
/**
 * This implementation uses recurrence relations.
 */
      template <class T>
        static T Pnm(int n1, int m1, const T &theta)
      {
        T P00(int(1));
// cf is needed for negative m.
        double cf;
        int n, m;
        if (n1>=0)
        {
          n=n1;
        }
        else
        {
          n=-n1-1;
        }
        m=(int)std::abs(m1);
        if (n==0 && m==0)
        {
          return P00;
        }
        if (m>n)
        {
          return T(int(0));
        }
// NOTE: in ginac and gmp maybe it's better to multiply at the end and generate
// cf on-the-fly, so that we can "mult_by_frac".
        if (m1>=0)
        {
          cf=1.;
        }
        else
        {
          cf=((double)cs_phase(m)*factorial(n-m))/factorial(n+m);
        }
        T retval(P00), oldPnm(int(0)), tmp1;
        int i;
        const std::complex<T> expitheta=complexp(theta);
        const T costheta=expitheta.real(), sintheta=expitheta.imag();
// Recursion to get from P_00 to P_mm.
        for (i=0;i<m;++i)
        {
          retval*=(-(2*i+1));
          retval*=sintheta;
        }
// Recursion to get from P_mm to P_nm (n>m).
        for (i=m;i<n;++i)
        {
// NOTE: in ginac and gmp we may use a "mult_by_frac" method to avoid losing
// arbitrary precision. Also in Pnm_vector below.
          oldPnm*=(-i-m)/(i+1.-m);
          tmp1=oldPnm;
          oldPnm=retval;
          retval*=(2*i+1)/(i+1.-m);
          retval*=costheta;
          retval+=tmp1;
        }
        retval*=cf;
        return retval;
      }
/// Non-normalized spherical harmonic.
      template <class T>
        static std::complex<T> Ynm(int n, int m, const T
        &theta, const T &phi)
      {
        T m_phi=phi;
        m_phi*=m;
        return (complexp(m_phi)*Pnm(n,m,theta));
      }
/// Vector of Pnm(cos(theta)) from P-n,n to Pnn.
      template <class T>
        static std::vector<T> Pnm_vector(int n1, const T &theta)
      {
        int n, i, j;
        if (n1>=0)
        {
          n=n1;
        }
        else
        {
          n=-n1-1;
        }
        std::vector<T> Pmm_vec(n+1);
        Pmm_vec[0]=T(1);
        if (n==0)
        {
          return Pmm_vec;
        }
        const std::complex<T> expitheta=complexp(theta);
        const T costheta=expitheta.real(), sintheta=expitheta.imag();
// Create vector of Pmm (0<=m<=n).
        for (i=1;i<n+1;++i)
        {
          Pmm_vec[i]=Pmm_vec[i-1];
// Here it is (i-1) instead of just i like above because in the sum the first
// term is 0, not 1
          Pmm_vec[i]*=(-(2*(i-1)+1));
          Pmm_vec[i]*=sintheta;
        }
        std::vector<T> retval(2*n+1);
        int m;
        for (i=n;i<2*n+1;++i)
        {
// Assign Pii.
          m=i-n;
          T oldPnm(int(0)), tmp1;
          retval[i]=Pmm_vec[m];
          T &ret=retval[i];
          for (j=m;j<n;++j)
          {
            oldPnm*=(-j-m)/(j+1.-m);
            tmp1=oldPnm;
            oldPnm=ret;
            ret*=(2*j+1)/(j+1.-m);
            ret*=costheta;
            ret+=tmp1;
          }
        }
// Make some space.
        Pmm_vec.clear();
// Take care of negative m.
        double cf;
        for (i=n-1;i>=0;--i)
        {
          m=i-n;
// -m instead of m here, because m<0.
          cf=((double)cs_phase(-m)*factorial(n+m))/factorial(n-m);
          retval[i]=retval[n-m];
          retval[i]*=cf;
        }
        return retval;
      }
/// Vector of spherical harmonics from Yn,-m to Yn,m.
      template <class T>
        static std::vector<std::complex<T> >
        Ynm_vector(int n, const T &theta, const T &phi)
      {
        std::vector<T> Pnm_vec=Pnm_vector(n,theta);
        std::vector<std::complex<T> > retval(2*n+1);
        int m;
        for (unsigned int i=0;i<retval.size();++i)
        {
          m=i-n;
          T m_phi=phi;
          m_phi*=m;
          retval[i]=complexp(m_phi)*Pnm_vec[i];
        }
        return retval;
      }
/// Vector of Bessel functions of the first kind from 0 to n.
// FIXME.
      template <class T>
        static std::vector<T> besselJ_vector(unsigned int n, const T &x)
      {
        std::vector<T> retval((size_t)(n+1));
        if (n==0)
        {
          retval[0]=besselJ(0,x);
          return retval;
        }
        else if (n==1)
        {
          retval[0]=besselJ(0,x);
          retval[1]=besselJ(1,x);
          return retval;
        }
        retval[0]=besselJ(0,x);
        retval[1]=besselJ(1,x);
        for (unsigned int j=2;j<=n;++j)
        {
          retval[j]=2*(j-1)*(retval[j-1]/x)-retval[j-2];
        }
        return retval;
      }
/// dnkm of Wigner's theorem.
/**
 * Input is the complex exponential of the angle.
 */
      template <class T>
        static T dnkm(int n, int k, int m, const std::complex<T> &input)
      {
        T retval(int(0)), tmp;
        double cf;
        const int t1=max(0,k-m), t2=min(n-m,n+k);
        const T real_input=input.real(), imag_input=input.imag();
        for (int t=t1;t<=t2;++t)
        {
          cf=(double)(cs_phase(t)*factorial(n-k)*factorial(n+m))/(factorial(t)*factorial(n+k-t)*
            factorial(n-m-t)*factorial(m-k+t));
// Assert >=0, we are going to cast to unsigned int so better safe than sorry.
          p_assert(2*n-m+k-2*t>=0);
          p_assert(m-k+2*t>=0);
// FIXME: cache the results here? There are lots of natural powers here.
          tmp=natural_pow((unsigned int)(2*n-m+k-2*t),real_input)*natural_pow((unsigned
            int)(m-k+2*t),imag_input);
          tmp*=cf;
          retval+=tmp;
        }
        return retval;
      }
/// Quadratic phase of integer n.
      static std::complex<double> quad_phase(int n)
      {
        if (is_odd(n))
        {
          if (is_odd((n-1)<<1))
          {
            return std::complex<double>(0.,-1.);
          }
          else
          {
            return std::complex<double>(0.,1.);
          }
        }
        else
        {
          if (is_odd(n<<1))
          {
            return std::complex<double>(-1.);
          }
          else
          {
            return std::complex<double>(1.);
          }
        }
      }
/// Dnkm of Wigern's theorem.
      template <class T>
        static std::complex<T>
        Dnkm(int n, int k, int m, const T &alpha,
        const std::complex<T> &complexhalfbeta, const T &gamma)
      {
        T dnkm_tmp(dnkm<T>(n,k,m,complexhalfbeta));
        return complexp(alpha*(-k)-gamma*(m))*(
          std::complex<T>(quad_phase(k-m))*=dnkm_tmp);
      }
/// Wigner rotation theorem.
/**
 * Returns Ynm'(theta,phi), i.e. Ynm(theta,phi) rotated under the three
 * Euler angles (alpha,beta,gamma).
 */
      template <class T>
        static std::complex<T> wig_rot(int n, int m,
        const T &alpha, const T &beta, const T &gamma, const T &theta, const T &phi)
      {
// NOTE: what does it give to multiply by .5 a series with non-zero linargs?
        std::complex<T> complexhalfbeta=complexp(beta*.5);
        std::vector<std::complex<T> > Ynm_vec=Ynm_vector(n,theta,phi);
        std::complex<T> retval(0.,0.);
        std::complex<T> tmp;
        for (int k=-n;k<=n;++k)
        {
          tmp=Ynm_vec[k+n];
          tmp*=Dnkm(n,k,m,alpha,complexhalfbeta,gamma);
          retval+=tmp;
        }
        return retval;
      }
/// Natural power.
      template <class T>
        static T natural_pow(unsigned int n, const T &x)
      {
        if (n==0)
        {
          return T(int(1));
        }
        T retval=x;
        for (unsigned int i=1;i<n;++i)
        {
          retval*=x;
        }
        return retval;
      }
/// Complex exponential, generic version.
      template <class T>
        static std::complex<T> complexp(const T &x)
      {
        return x.complexp();
      }
/// Cosine, generic version.
      template <class T>
        static T cosine(const T &x)
      {
        return x.cosine();
      }
/// Sine, generic version.
      template <class T>
        static T sine(const T &x)
      {
        return x.sine();
      }
/**
 * Return absolute value of difference from nearest integer.
 * @param x input value.
 */
      static double delta_nearbyint(const double &x)
      {
        return std::abs(x-nearbyint(x));
      }
    private:
/// Private ctor.
      math()
        {}
// Low-level factorial function.
      template <class T>
        static T ll_factorial(int n)
      {
        if (n<0)
        {
          std::cout << "FATAL: factorial of a negative number." << std::cout;
          std::exit(1);
        }
        T k=1;
        int i=1;
        while (i<n)
        {
          ++i;
          k*=i;
        }
        return k;
      }
  };

/// Find norm, specialization for double.
  template <>
    inline double math::norm<double>(const double &x)
  {
    return std::abs(x);
  }

/// Complex exponential, specialization for double.
  template <>
    inline std::complex<double> math::complexp<double>(const double &x)
  {
    return std::polar(1.,x);
  }

/// Cosine, specialization for double.
  template <>
    inline double math::cosine<double>(const double &x)
  {
    return std::cos(x);
  }

/// Sine, specialization for double.
  template <>
    inline double math::sine<double>(const double &x)
  {
    return std::sin(x);
  }
}
#endif
