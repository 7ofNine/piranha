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

#ifndef PIRANHA_BASE_PSERIES_PROBE_MP_H
#define PIRANHA_BASE_PSERIES_PROBE_MP_H

/** @file base_pseries_probe_mp.h
    \brief Meta-programming for series probing.
*/

#include "../../piranha_tbb.h" // For parallel evaluation.

namespace piranha
{
/// Base class for the evaluation of a series over a timespan.
/**
 * Upon construction this class will simply perform some sanity checks on the interval evaluation parameters,
 * set a suitability flag and then return.
 */
  template <class Series>
    class base_series_interval_evaluator
  {
      typedef typename Series::eval_type eval_type;
    public:
/// Constructor from interval parameters.
/**
 * If the parameters are sane, is_suitable will be set to true and the retval vector will be appropriately resized,
 * otherwise is_suitable will be set to false.
 */
      base_series_interval_evaluator(const Series &s, const double &t0, const double &t1,
        const int &n, std::vector<eval_type> &retval):is_suitable(true),m_series(s),m_t0(t0),
        m_t1(t1),m_n(n),m_retval(retval)
      {
        preliminary_checks();
        if (!is_suitable)
        {
          return;
        }
        m_size = (size_t)n;
        m_retval.clear();
        m_retval.resize(m_size);
      }
      void serial_evaluation() const
      {
        if (!is_suitable)
        {
          return;
        }
        double t=m_t0;
        for (size_t i=0;i < m_size;++i)
        {
          m_retval[i]=m_series.t_eval(t);
          t+=m_step;
        }
      }
/// Return is_suitable flag.
      bool status() const {return is_suitable;}
    private:
      base_series_interval_evaluator() {}
      void preliminary_checks()
      {
        if (m_n <= 0)
        {
          std::cout << "Please insert a strictly positive value for the number of steps in interval series evaluation."
            << std::endl;
          is_suitable=false;
          return;
        }
        m_step = (m_t1-m_t0)/(double)m_n;
// Check that step is not null and that signs of interval and step are consistent.
        if (m_step == 0 or (m_t1-m_t0) * m_step < 0)
        {
          std::cout << "Error: problem in step size in interval series evaluation." << std::endl;
          is_suitable=false;
          return;
        }
      }
    public:
/// Step of the interval.
      double                  m_step;
/// Suitability flag.
      bool                    is_suitable;
/// Size of the interval.
      size_t                  m_size;
/// Const reference to series.
      const Series            &m_series;
/// Const reference to starting point of time interval.
      const double            &m_t0;
/// Const reference to ending point of time interval.
      const double            &m_t1;
/// Requested size of the interval.
      int                     m_n;
/// Reference to return value.
      std::vector<eval_type>  &m_retval;
  };

/// Serial evaluation of a series over a timespan.
  template <bool Parallel, class Series>
    class series_interval_evaluator:public base_series_interval_evaluator<Series>
  {
      typedef typename Series::eval_type eval_type;
      typedef base_series_interval_evaluator<Series> ancestor;
    public:
      series_interval_evaluator(const Series &s, const double &t0, const double &t1, const int &n,
        std::vector<eval_type> &retval):ancestor::base_series_interval_evaluator(s,t0,t1,n,retval)
      {
        ancestor::serial_evaluation();
      }
    private:
      series_interval_evaluator() {}
  };

#ifdef _PIRANHA_TBB
/// Parallel evaluation of a series over a timespan.
  template <class Series>
    class series_interval_evaluator<true,Series>:public base_series_interval_evaluator<Series>
  {
      typedef typename Series::eval_type eval_type;
      typedef base_series_interval_evaluator<Series> ancestor;
      typedef series_interval_evaluator<true,Series> sie_type;
      class parallel_series_evaluation
      {
        public:
          parallel_series_evaluation(const sie_type &sie):m_sie(sie) {}
          void operator()(const tbb::blocked_range<size_t> &r) const
          {
            double t=m_sie.m_t0;
            for(size_t i=r.begin();i != r.end();++i)
            {
              m_sie.m_retval[i]=m_sie.m_series.t_eval(t);
              t+=m_sie.m_step;
            }
          }
        private:
          parallel_series_evaluation() {}
          const sie_type  &m_sie;
      };
    public:
      series_interval_evaluator(const Series &s, const double &t0, const double &t1, const int &n,
        std::vector<eval_type> &retval):ancestor::base_series_interval_evaluator(s,t0,t1,n,retval)
      {
        if (!ancestor::is_suitable)
        {
          return;
        }
// Parallel version.
        tbb::parallel_for(tbb::blocked_range<size_t>(0,ancestor::m_size,100),parallel_series_evaluation(*this));
      }
    private:
      series_interval_evaluator() {}
  };
#endif
}

#endif
