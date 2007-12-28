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

#ifndef PIRANHA_PIRANHA_TBB_H
#define PIRANHA_PIRANHA_TBB_H

#ifdef _PIRANHA_TBB

#include <tbb/blocked_range.h>
#include <tbb/task_scheduler_init.h>
#include <tbb/parallel_for.h>
#include <vector>

namespace piranha
{
  static const tbb::task_scheduler_init tbb_init;

/// Parallel evaluation of a series over a timespan.
  template <class Series>
    class parallel_series_evaluation
  {
      typedef typename Series::eval_type eval_type;
    public:
      parallel_series_evaluation(std::vector<eval_type> &retval, const Series &series,
        const double &step, const double &t0):m_retval(retval),m_series(series),m_step(step),m_t0(t0) {}
      void operator()(const tbb::blocked_range<size_t> &r) const
      {
        std::vector<eval_type> *retval=&m_retval;
        Series const *series=&m_series;
        const double step = m_step;
        double t=m_t0;
        for(size_t i=r.begin();i != r.end();++i)
        {
          (*retval)[i]=series->t_eval(t);
          t+=step;
        }
      }
    private:
      std::vector<eval_type>  &m_retval;
      const Series            &m_series;
      const double            &m_step;
      const double            &m_t0;
  };
}

#endif

#endif
