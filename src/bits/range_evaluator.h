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

#ifndef PIRANHA_RANGE_EVALUATOR_H
#define PIRANHA_RANGE_EVALUATOR_H

#include "compile_switches.h"                     // For parallel mode.
#include "piranha_tbb.h"                          // For parallel evaluation.

namespace piranha
{
  /// Base class for the evaluation over a range.
  /**
   * Upon construction this class will simply perform some sanity checks on the interval evaluation parameters
   * and then return. In order to be evaluated, a class must conform to the
   * piranha::concepts::evaluatable concept.
   */
  template <class Evaluatable>
    class base_range_evaluator
  {
    protected:
      typedef typename Evaluatable::eval_type eval_type;
    public:
      /// Constructor from interval parameters.
      /**
       * If the parameters are sane the retval vector will be appropriately resized, otherwise
       * size will be zero.
       */
      base_range_evaluator(const Evaluatable &e, const double &t0, const double &t1, const int &n):
      m_evaluated(e),m_t0(t0),m_t1(t1),m_n(n),m_size(0),m_step(0),m_retval(),m_time()
      {
        preliminary_checks();
        m_size = (size_t)n;
        m_retval.resize(m_size);
        m_time.resize(m_size);
      }
      /// Return const reference to return values vector.
      const std::vector<eval_type> &values() const {return m_retval;}
      /// Return const reference to interval vector.
      const std::vector<double> &times() const {return m_time;}
    protected:
      /// Serial evaluation.
      void serial_evaluation()
      {
        double t=m_t0;
        for (size_t i=0;i < m_size;++i)
        {
          m_retval[i]=m_evaluated.t_eval(t);
          m_time[i]=t;
          t+=m_step;
        }
      }
    private:
      base_range_evaluator() {}
      void preliminary_checks()
      {
        if (m_n <= 0)
        {
          std::cout << "Please insert a strictly positive value for the number of steps in range evaluation."
            << std::endl;
          return;
        }
        m_step = (m_t1-m_t0)/(double)m_n;
        // Check that step is not null and that signs of interval and step are consistent.
        if (m_step == 0 or (m_t1-m_t0) * m_step < 0)
        {
          std::cout << "Error: problem in step size in range evaluation." << std::endl;
          return;
        }
      }
    public:
      // These are all references because they are used only in contruction.
      /// Const reference to series.
      const Evaluatable       &m_evaluated;
      /// Const reference to starting point of time interval.
      const double            &m_t0;
      /// Const reference to ending point of time interval.
      const double            &m_t1;
      /// Requested size of the interval.
      const int               &m_n;
      /// Size of the interval.
      size_t                  m_size;
      /// Step of the interval.
      double                  m_step;
      /// Return values.
      std::vector<eval_type>  m_retval;
      /// Time array.
      std::vector<double>     m_time;
  };

  /// Serial range evaluation.
  template <class Evaluatable, bool Parallel = compile_switches::use_tbb>
    class range_evaluator:public base_range_evaluator<Evaluatable>
  {
    typedef base_range_evaluator<Evaluatable> ancestor;
    typedef typename ancestor::eval_type eval_type;
    public:
      range_evaluator(const Evaluatable &e, const double &t0, const double &t1, const int &n):
      ancestor::base_range_evaluator(e,t0,t1,n) {ancestor::serial_evaluation();}
    private:
      range_evaluator() {}
  };

#ifdef _PIRANHA_TBB
  /// Parallel range evaluation.
  template <class Evaluatable>
    class range_evaluator<Evaluatable,true>:public base_range_evaluator<Evaluatable>
  {
    typedef base_range_evaluator<Evaluatable> ancestor;
    typedef typename ancestor::eval_type eval_type;
    typedef range_evaluator<Evaluatable,true> re_type;
    class parallel_range_evaluation
    {
      public:
        parallel_range_evaluation(re_type &re):m_re(re) {}
        void operator()(const tbb::blocked_range<size_t> &r) const
        {
          for(size_t i=r.begin();i != r.end();++i)
          {
            const double t = m_re.m_t0+i*m_re.m_step;
            m_re.m_retval[i]=m_re.m_evaluated.t_eval(t);
            m_re.m_time[i]=t;
          }
        }
      private:
        parallel_range_evaluation() {}
        re_type &m_re;
    };
    public:
      range_evaluator(const Evaluatable &e, const double &t0, const double &t1, const int &n):
      ancestor::base_range_evaluator(e,t0,t1,n)
      {
                                                  // HARDCODED!
        const size_t grain_size = ancestor::m_size/100;
        switch (grain_size == 0)
        {
          case true:
            ancestor::serial_evaluation();
            break;
          case false:
            // Parallel version.
            tbb::parallel_for(tbb::blocked_range<size_t>(0,ancestor::m_size,grain_size),parallel_range_evaluation(*this));
        }
      }
    private:
      range_evaluator() {}
  };
#endif
}
#endif
