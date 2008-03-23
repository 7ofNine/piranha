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

#ifndef PIRANHA_SERIES_MULTIPLIER_MP_H
#define PIRANHA_SERIES_MULTIPLIER_MP_H

#include <boost/tuple/tuple.hpp>

namespace piranha
{
  struct base_insert_multiplication_result
  {
    protected:
      template <class SingleRes, class MultSet, class ArgsTuple, class Truncator>
        static void insert_single_res(const SingleRes &res, MultSet &mult_set, const ArgsTuple &args_tuple,
        const Truncator &trunc)
      {
        const typename MultSet::const_iterator it = mult_set.find(res);
        switch (trunc.accept_term(res))
        {
          case true:
            switch (it == mult_set.end())
            {
              // Not duplicate, insert it.
              case true:
                mult_set.insert(res);
                break;
              // Duplicate, modify existing.
              case false:
                it->m_cf.add(res.m_cf,args_tuple);
            }
            break;
          case false:
            ;
        }
      }
  };

  // This gets specialised here and below to handle both the case in which the multiplication result
  // is a single entity, and that in which multiplication result is a cons list of return values.
  // This is needed to handle the case of Poisson series, in which we have two terms resulting from each
  // term-by-term multiplication.
  template <class SingleRes>
    struct insert_multiplication_result:public base_insert_multiplication_result
  {
    template <class MultSet, class ArgsTuple, class Truncator>
      static void run(const SingleRes &res, MultSet &mult_set, const ArgsTuple &args_tuple, const Truncator &trunc)
    {
      insert_single_res(res,mult_set,args_tuple,trunc);
    }
  };

  template <class Head, class Tail>
    struct insert_multiplication_result<boost::tuples::cons<Head,Tail> >:public base_insert_multiplication_result
  {
    template <class MultSet, class ArgsTuple, class Truncator>
      static void run(const boost::tuples::cons<Head,Tail> &mult_res, MultSet &mult_set, const ArgsTuple &args_tuple,
      const Truncator &trunc)
    {
      insert_single_res(mult_res.get_head(),mult_set,args_tuple,trunc);
      insert_multiplication_result<Tail>::run(mult_res.get_tail(),mult_set,args_tuple,trunc);
    }
  };

  template <class Head>
    struct insert_multiplication_result<boost::tuples::cons<Head,boost::tuples::null_type> >:
    public base_insert_multiplication_result
  {
    template <class MultSet, class ArgsTuple, class Truncator>
      static void run(const boost::tuples::cons<Head,boost::tuples::null_type> &mult_res, MultSet &mult_set,
      const ArgsTuple &args_tuple, const Truncator &trunc)
    {
      insert_single_res(mult_res.get_head(),mult_set,args_tuple,trunc);
    }
  };
}

#endif
