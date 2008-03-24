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

#ifndef PIRANHA_FOURIER_SERIES_MULTIPLIER_H
#define PIRANHA_FOURIER_SERIES_MULTIPLIER_H

#include "../base_classes/plain_series_multiplier.h"

namespace piranha
{
  /// Series multiplier specifically tuned for Fourier series.
  /**
   * This multiplier internally will used coded arithmetics if possible, otherwise it will operate just
   * like piranha::plain_series_multiplier.
   */
  template <class Series1, class Series2, class ArgsTuple, template <class, class, class> class Truncator>
    class fourier_series_multiplier:public plain_series_multiplier<Series1,Series2,ArgsTuple,Truncator>
  {
      typedef plain_series_multiplier<Series1,Series2,ArgsTuple,Truncator> ancestor;
      typedef typename ancestor::truncator_type truncator_type;
      typedef typename ancestor::term_type term_type;
      typedef typename ancestor::cf1 cf1;
      typedef typename ancestor::cf2 cf2;
      typedef typename ancestor::key key;
      typedef typename ancestor::sm_cf1 sm_cf1;
      typedef typename ancestor::sm_cf2 sm_cf2;
      typedef typename ancestor::sm_key sm_key;
      typedef typename ancestor::mult_set mult_set;
    public:
      fourier_series_multiplier(const Series1 &s1, const Series2 &s2, Series1 &retval, const ArgsTuple &args_tuple):
        ancestor::plain_series_multiplier(s1,s2,retval,args_tuple)
      {}
      /// Perform multiplication and place the result into m_retval.
      void perform_multiplication()
      {
        ancestor::adjust_sizes();
        ancestor::cache_series_terms(ancestor::m_s1,ancestor::m_cfs1,ancestor::m_keys1);
        ancestor::cache_series_terms(ancestor::m_s2,ancestor::m_cfs2,ancestor::m_keys2);
        ancestor::plain_multiplication();
        ancestor::plain_insert_result_into_retval();
      }
  };
}

#endif
