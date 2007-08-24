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

#ifndef PIRANHA_PRECTEST_H
#define PIRANHA_PRECTEST_H

#include <fstream>
#include <iostream>
#include <valarray>

#include "common_typedefs.h"                      // tc_add_ps_to_arg.
#include "math.h"
#include "phase_list.h"                           // tc_insert_phase.
#include "stream_manager.h"                       // Gnuplot save.
#include "utils.h"                                // Check filename dir.

namespace piranha
{
/// Term-by-term comparison of two series.
/**
 * The norms of the coefficients are compared on a term-by-term basis, and the absolute values of
 * the differences are stored internally. If an analogous term cannot be found a negative value
 * is produced, which is made proportional to the largest difference detected. The differences can
 * be relative (taking the first series' coefficients as base) or absolute.
 *
 * The results of the comparison can the be retrieved with the sc::diffs method, or saved to
 * text files to be used by gnuplot to display graphically the comparison.
 */
  template <class T>
    class sc
  {
    public:
      sc(const T &, const T &, bool relative=true);
/// Return differences.
      const double &diffs(const size_t &n) const
      {
        return diffs_[n];
      }
      void gnuplot_save(const std::string &) const;
/// Get size of diffs vector.
      size_t size() const
      {
        return diffs_.size();
      }
/// Get type of differences (absolute or relative).
// FIXME: use enum?
      bool is_relative() const
      {
        return is_relative_;
      }
    private:
      sc()
      {
        std::cout << "PRIVATE CTOR!" << std::endl;
        std::exit(1);
      }
// Data members
      std::valarray<double>   diffs_;
      bool                    is_relative_;
  };

/// Constructor from series.
  template <class T>
    inline sc<T>::sc(const T &ps1, const T &ps2, bool relative):
  diffs_(ps1.length()),is_relative_(relative)
  {
    if (ps1.length()==0)
    {
      return;
    }
    typename T::it_h_index it_h;
    const typename T::it_s_index it_f=ps1.g_s_index().end();
    size_t i=0;
    double max_diff=0.;
    for (typename T::it_s_index it=ps1.g_s_index().begin();it!=it_f;++it)
    {
      it_h=ps2.g_h_index().find(boost::make_tuple(it->g_flavour(),
        *it->g_trig()));
      if (it_h==ps2.g_h_index().end())
      {
        diffs_[i]=-.1;
      }
      else
      {
        if (relative)
        {
          diffs_[i]=std::abs(1.-std::abs(it_h->g_cf()->norm(ps2.cf_s_vec())/it->g_cf()->norm(ps1.cf_s_vec())));
        }
        else
        {
          diffs_[i]=std::abs(it_h->g_cf()->norm(ps2.cf_s_vec())-it->g_cf()->norm(ps1.cf_s_vec()));
        }
        if (diffs_[i]>max_diff)
        {
          max_diff=diffs_[i];
        }
      }
      ++i;
    }
    for (i=0;i<diffs_.size();++i)
    {
      if (diffs_[i]<0)
      {
        diffs_[i]=-max_diff/10.;
      }
    }
  }

/// Save results of the comparison in gnuplot format.
  template <class T>
    inline void sc<T>::gnuplot_save(const std::string &filename) const
  {
    const std::string plot_file=filename+".plt", data_file=filename+".dat";
    std::ofstream outf_plot(plot_file.c_str(),std::ios_base::trunc);
    std::ofstream outf_data(data_file.c_str(),std::ios_base::trunc);
    if (outf_plot.fail() || outf_data.fail())
    {
      std::cout << "Error saving to file " << filename << std::endl;
      outf_plot.close();
      outf_data.close();
      return;
    }
    stream_manager::setup_print(outf_data);
    for (size_t i=0;i<diffs_.size();++i)
    {
      outf_data << diffs_[i] << std::endl;
    }
    outf_plot << "set xrange [0:" << diffs_.size() << "]" << std::endl;
    outf_plot << "plot '" << data_file <<  "' with impulses";
    outf_plot << std::endl << "reset";
    outf_plot.close();
    outf_data.close();
  }

/// Series comparison in the time domain
/** Two series are evaluated over a time span and statistics on
 * the evaluations are computed, like maximum error and rms of error (sigma). Beside comparing two
 * series it is also possible to compare a series against an operation: for instance we can compare
 * the result of a multiplication against the "correct" result, i.e. the product of the evaluation
 * of the two factors of the multiplication over the given time span. Comparisons for most operations
 * are provided.
 */
  template <class T>
    class base_tc
  {
    public:
/// Alias for evaluation type.
      typedef typename T::eval_type eval_type;
/// Get number of time comparisons performed.
      size_t size() const
      {
        return time_.size();
      }
/// Return time at comparison n.
      const double &time(const size_t &n) const
      {
        return time_[n];
      }
/// Return value of resulting series at comparison n.
      const eval_type &hs(const size_t &n) const
      {
        return hs_[n];
      }
/// Return exact value at comparison n.
      const eval_type &hs_computed(const size_t &n) const
      {
        return hs_computed_[n];
      }
/// Return error at comparison n.
      double error(const size_t &n) const
      {
        return std::abs(hs(n)-hs_computed(n));
      }
/// Get rms.
      const double &sigma() const
      {
        return sigma_;
      }
/// Get max error.
      const double &max_error() const
      {
        return max_error_;
      }
      void gnuplot_save(const std::string &) const;
    protected:
/// Calculate diff at time t.
/**
 * This function will have to be re-implemented in inherited class to perform the desired
 * comparisons.
 */
      virtual eval_type eval_hs_computed(const double &t) const = 0;
      void init(const double &t1, const double &t2, const size_t &ntot,
        const T &benchmarked)
      {
        t1_=t1;
        t2_=t2;
        ntot_=ntot;
        time_.resize(ntot);
        hs_.resize(ntot);
        hs_computed_.resize(ntot);
        benchmarked_=&benchmarked;
        build_tc();
      }
// Ctor & Dtor
      base_tc():benchmarked_(0)
        {}
      ;
      virtual ~base_tc()
        {}
      ;
    private:
      eval_type eval_hs(const double &t) const
      {
        return benchmarked_->t_eval(t);
      }
      void build_tc();
      void calc_stats();
      void plot_data(std::ostream &) const;
// Data members
      double                      t1_;
      double                      t2_;
      size_t                      ntot_;
      std::valarray<double>       time_;
      std::valarray<eval_type>    hs_;
      std::valarray<eval_type>    hs_computed_;
      double                      sigma_;
      double                      max_error_;
      const T                     *benchmarked_;
  };

  template <class T>
    inline void base_tc<T>::build_tc()
  {
    const double step_size=(t2_-t1_)/ntot_;
    double t=t1_;
    for (size_t i=0;i<ntot_;++i)
    {
      time_[i]=t;
      hs_[i]=eval_hs(t);
      hs_computed_[i]=eval_hs_computed(t);
      t+=step_size;
    }
    std::cout << "Computing statistics..." << std::endl;
    calc_stats();
    std::cout << "Done." << std::endl;
  }

  template <class T>
    inline void base_tc<T>::calc_stats()
  {
    double tmp=0., max=0., candidate;
    if (hs_.size()==0)
    {
      sigma_=tmp;
      max_error_=max;
      return;
    }
    size_t i;
    for (i=0;i<hs_.size();++i)
    {
      candidate=std::abs(hs_[i]-hs_computed_[i]);
      if (candidate>max)
      {
        max=candidate;
      }
      tmp+=math::natural_pow(2,candidate);
    }
    sigma_=std::sqrt(tmp/(i+1.));
    max_error_=max;
  }

/// Save results in gnuplot format.
  template <class T>
    inline void base_tc<T>::gnuplot_save(const std::string &filename) const
  {
    const std::string plot_file=filename+".plt", data_file=filename+".dat";
    std::ofstream outf_plot(plot_file.c_str(),std::ios_base::trunc);
    std::ofstream outf_data(data_file.c_str(),std::ios_base::trunc);
    if (outf_plot.fail() || outf_data.fail())
    {
      std::cout << "Error saving to file " << filename << std::endl;
      outf_plot.close();
      outf_data.close();
      return;
    }
    stream_manager::setup_print(outf_data);
// Setup plot commands
    outf_plot << "set logscale y" << std::endl;
    outf_plot << "plot '" << data_file << "' using 1:2 title 'diff', '" <<
      data_file << "' using 1:3 title 'exact'";
    outf_plot << std::endl << "reset";
    outf_plot.close();
// Write to data file
    plot_data(outf_data);
    outf_data.close();
  }

// Print plot data to file

  template <class T>
    inline void base_tc<T>::plot_data(std::ostream &os) const
  {
    for (size_t i=0;i<hs_.size();++i)
    {
      os << time_[i] << '\t' << std::abs(hs_[i]-hs_computed_[i]) << '\t' <<
        std::abs(hs_computed_[i]) << std::endl;
    }
  }

// Comparisons for math operations

  template <class T>
    class tc_equal:public base_tc<T>
  {
    public:
// b_type stands for "benchmarked type"
      typedef T b_type;
      typedef typename base_tc<T>::eval_type eval_type;
      tc_equal(const T &b, const double &t1, const double &t2, const size_t &ntot,
        const T &a):a_(&a)
      {
        base_tc<T>::init(t1,t2,ntot,b);
      }
    private:
      const T     *a_;
      virtual eval_type eval_hs_computed(const double &t) const
      {
        return a_->t_eval(t);
      }
  };

  template <class T>
    class tc_mult:public base_tc<T>
  {
    public:
      typedef T b_type;
      typedef typename base_tc<T>::eval_type eval_type;
      tc_mult(const T &b, const double &t1, const double &t2, const size_t &ntot,
        const T &x, const T &y):x_(&x),y_(&y)
      {
        base_tc<T>::init(t1,t2,ntot,b);
      }
    private:
      const T     *x_;
      const T     *y_;
      virtual eval_type eval_hs_computed(const double &t) const
      {
        return x_->t_eval(t)*y_->t_eval(t);
      }
  };

  template <class T>
    class tc_complexp:public base_tc<std::complex<T> >
  {
    public:
      typedef std::complex<T> b_type;
      typedef typename base_tc<b_type>::eval_type eval_type;
      tc_complexp(const b_type &b, const double &t1, const double &t2,
        const size_t &ntot, const T &a):a_(&a)
      {
        base_tc<b_type>::init(t1,t2,ntot,b);
      }
    private:
      const T     *a_;
      virtual eval_type eval_hs_computed(const double &t) const
      {
        return math::complexp(a_->t_eval(t));
      }
  };

  template <class T>
    class tc_cosine:public base_tc<T>
  {
    public:
      typedef T b_type;
      typedef typename base_tc<T>::eval_type eval_type;
      tc_cosine(const T &b, const double &t1, const double &t2,
        const size_t &ntot, const T &a):a_(&a)
      {
        base_tc<T>::init(t1,t2,ntot,b);
      }
    private:
      const T     *a_;
      virtual eval_type eval_hs_computed(const double &t) const
      {
        return std::cos(a_->t_eval(t));
      }
  };

  template <class T>
    class tc_sine:public base_tc<T>
  {
    public:
      typedef T b_type;
      typedef typename base_tc<T>::eval_type eval_type;
      tc_sine(const T &b, const double &t1, const double &t2,
        const size_t &ntot, const T &a):a_(&a)
      {
        base_tc<T>::init(t1,t2,ntot,b);
      }
    private:
      const T     *a_;
      virtual eval_type eval_hs_computed(const double &t) const
      {
        return std::sin(a_->t_eval(t));
      }
  };

  template <class T>
    class tc_Pnm:public base_tc<T>
  {
    public:
      typedef T b_type;
      typedef typename base_tc<T>::eval_type eval_type;
      tc_Pnm(const T &b, const double &t1, const double &t2,
        const size_t &ntot, int n, int m, const T &a):a_(&a),n_(n),m_(m)
      {
        base_tc<T>::init(t1,t2,ntot,b);
      }
    private:
      const T     *a_;
      int         n_;
      int         m_;
      virtual eval_type eval_hs_computed(const double &t) const
      {
        return math::Pnm(n_,m_,a_->t_eval(t));
      }
  };

  template <class T>
    class tc_Ynm:public base_tc<std::complex<T> >
  {
    public:
      typedef std::complex<T> b_type;
      typedef typename base_tc<b_type>::eval_type eval_type;
      tc_Ynm(const b_type &b, const double &t1, const double &t2,
        const size_t &ntot, int n, int m, const T &theta, const T
        &phi):theta_(&theta),phi_(&phi),
        n_(n),m_(m)
      {
        base_tc<b_type>::init(t1,t2,ntot,b);
      }
    private:
      const T     *theta_;
      const T     *phi_;
      int         n_;
      int         m_;
      virtual eval_type eval_hs_computed(const double &t) const
      {
        return math::Ynm(n_,m_,theta_->t_eval(t),phi_->t_eval(t));
      }
  };

  template <class T>
    class tc_wig_rot:public base_tc<std::complex<T> >
  {
    public:
      typedef std::complex<T> b_type;
      typedef typename base_tc<b_type>::eval_type eval_type;
      tc_wig_rot(const b_type &b, const double &t1,
        const double &t2,const size_t &ntot, int n, int m, const T &alpha,
        const T &beta, const T &gamma,const T &theta, const T &phi):n_(n),
        m_(m),alpha_(&alpha),beta_(&beta),gamma_(&gamma),theta_(&theta),phi_(&phi)
      {
        base_tc<b_type>::init(t1,t2,ntot,b);
      }
    private:
      int         n_;
      int         m_;
      const T     *alpha_;
      const T     *beta_;
      const T     *gamma_;
      const T     *theta_;
      const T     *phi_;
      virtual eval_type eval_hs_computed(const double &t) const
      {
        return math::wig_rot(n_,m_,alpha_->t_eval(t),beta_->t_eval(t),gamma_->t_eval(t),
          theta_->t_eval(t),phi_->t_eval(t));
      }
  };

  template <class T>
    class tc_pow:public base_tc<T>
  {
    public:
      typedef T b_type;
      typedef typename base_tc<T>::eval_type eval_type;
      tc_pow(const T &b, const double &t1, const double &t2,
        const size_t &ntot, const double &power_, const T &a):a_(&a),power(power_)
      {
        base_tc<T>::init(t1,t2,ntot,b);
      }
    private:
      const T     *a_;
      double      power;
      virtual eval_type eval_hs_computed(const double &t) const
      {
        return std::pow(a_->t_eval(t),power);
      }
  };

// FIXME: here and below it would be better to use a function defined in series, instead of
// duplicating code.
/// Time comparison for addition of a series to an argument.
/**
 * @see piranha::base_pseries::add_ps_to_arg. FIXME!
 */
  template <class T>
    class tc_add_ps_to_arg:public base_tc<T>
  {
    public:
      typedef T b_type;
      typedef typename base_tc<T>::eval_type eval_type;
      tc_add_ps_to_arg(const T &b, const double &t1, const double &t2,
        const size_t &ntot, std::string name, const T &a, const T &orig):a_(&a),orig_(&orig)
      {
        size_t idx=orig_->trig_index(name);
        sym_index_=idx;
        base_tc<T>::init(t1,t2,ntot,b);
      }
    private:
      const T     *a_;
      const T     *orig_;
      trig_size_t sym_index_;
      virtual eval_type eval_hs_computed(const double &t) const
      {
        typedef typename T::r_it_s_index r_it_s_index;
        eval_type retval(0.);
        const size_t w=orig_->trig_width();
        int multiplier;
        double tmp;
        eval_type c_eval;
// Terms - start from the smallest, so that we keep good precision in the summation
        const r_it_s_index it_f=orig_->g_s_index().rend();
        for (r_it_s_index it=orig_->g_s_index().rbegin();it!=it_f;++it)
        {
// Check that we are not going outside the boundaries.
          if (sym_index_ < orig_->trig_width())
          {
            multiplier=it->g_trig()->multiplier(sym_index_);
          }
          else
          {
            multiplier=0;
          }
          tmp=it->g_trig()->t_eval(t,orig_->trig_s_vec());
          c_eval=it->g_cf()->t_eval(t,orig_->cf_s_vec());
          switch (it->g_flavour())
          {
            case true:
              retval+=c_eval*std::cos(tmp)*std::cos(multiplier*a_->t_eval(t))-
                c_eval*std::sin(tmp)*std::sin(multiplier*a_->t_eval(t));
              break;
            case false:
              retval+=c_eval*std::sin(tmp)*std::cos(multiplier*a_->t_eval(t))+
                c_eval*std::cos(tmp)*std::sin(multiplier*a_->t_eval(t));
          }
        }
// Linear arguments
        for (size_t j=0;j<w;++j)
        {
          retval+=orig_->lin_args()[j]*orig_->trig_s_vec()[j]->t_eval(t);
        }
        return retval;
      }
  };

/// Time comparison for phase insertion.
/**
 * @see piranha::base_pseries::insert_phases.
 */
  template <class T>
    class tc_insert_phases:public base_tc<T>
  {
    public:
      typedef T b_type;
      typedef typename base_tc<T>::eval_type eval_type;
      tc_insert_phases(const T &b, const double &t1, const double &t2,
        const size_t &ntot, const phase_list &pl, const T &a):a_(&a),pl_(&pl)
      {
        base_tc<T>::init(t1,t2,ntot,b);
      }
    private:
      const T             *a_;
      const phase_list    *pl_;
      virtual eval_type eval_hs_computed(const double &t) const
      {
        typedef typename T::iterator iterator;
        eval_type retval(0.);
        const size_t w=a_->trig_width();
        phase_list::const_iterator it2=pl_->begin();
        double tmp_phase;
        double tmp;
        eval_type c_eval;
// Terms - start from the smallest, so that we keep good precision in the summation.
        const iterator it_f=a_->end();
        for (iterator it=a_->begin();it!=it_f;++it)
        {
          if (it2!=pl_->end())
          {
            switch (pl_->operation())
            {
              case phase_list::add
                :
              tmp_phase=*it2;
              break;
              default:
                tmp_phase=*it2-it->g_trig()->phase(a_->trig_s_vec());
            }
            tmp=it->g_trig()->t_eval(t,a_->trig_s_vec());
            c_eval=it->g_cf()->t_eval(t,a_->cf_s_vec());
            switch (it->g_flavour())
            {
              case true:
                retval+=c_eval*std::cos(tmp)*std::cos(tmp_phase)-
                  c_eval*std::sin(tmp)*std::sin(tmp_phase);
                break;
              case false:
                retval+=c_eval*std::sin(tmp)*std::cos(tmp_phase)+
                  c_eval*std::cos(tmp)*std::sin(tmp_phase);
            }
            ++it2;
          }
          else
          {
            retval+=it->t_eval(t,a_->cf_s_vec(),a_->trig_s_vec());
          }
        }
// Linear arguments.
        for (size_t j=0;j<w;++j)
        {
          retval+=a_->lin_args()[j]*a_->trig_s_vec()[j]->t_eval(t);
        }
        return retval;
      }
  };
}
#endif
