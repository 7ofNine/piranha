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

#ifndef PIRANHA_IPOLY_H
#define PIRANHA_IPOLY_H

#include "p_assert.h"

#include <boost/integer_traits.hpp>
#include <boost/static_assert.hpp>
#include <boost/type_traits/is_integral.hpp>
#include <boost/type_traits/is_pod.hpp>
#include <vector>

namespace piranha
{
// TODO: move to separate file with other cachers?
/// Natural power cacher for integral type.
  template <class T>
    class integral_npow_cache
  {
      BOOST_STATIC_ASSERT(boost::is_integral<T>::value);
      BOOST_STATIC_ASSERT(boost::is_pod<T>::value);
      BOOST_STATIC_ASSERT(!(boost::integer_traits<T>::is_signed));
      typedef std::vector<std::vector<T> > container;
    public:
      static const T &request(const int &n, const T &arg)
      {
        p_assert(n >= 0);
        while (arg >= private_cache_.size())
        {
// Add a row to the matrix.
          private_cache_.push_back(std::vector<T>());
// Add the first element to the row.
          private_cache_.back().push_back(1);
        }
        while ((size_t)n >= private_cache_[arg].size())
        {
//std::cout << "before: " << private_cache_[arg].back() << '\n';
          private_cache_[arg].push_back(arg*private_cache_[arg].back());
//std::cout << "now: " << private_cache_[arg].back() << '\n';
        }
        return private_cache_[arg][n];
      }
    private:
      static container    private_cache_;
  };

  template <class T>
    typename integral_npow_cache<T>::container integral_npow_cache<T>::private_cache_;

/// Indexed monomial.
  template <class Cf, class Index>
    class imonomial
  {
      BOOST_STATIC_ASSERT(boost::is_integral<Index>::value);
      BOOST_STATIC_ASSERT(boost::is_pod<Index>::value);
      BOOST_STATIC_ASSERT(!(boost::integer_traits<Index>::is_signed));
    public:
      imonomial()
        {}
      imonomial(const Cf &value, const Index &i):private_index_(i),private_cf_(value)
        {}
      ~imonomial()
        {}
      static const Index &g_max_index()
      {
        return private_max_index_;
      }
    private:
      Index               private_index_;
      Cf                  private_cf_;
      static const Index  private_max_index_ = ((((((Index)1)<<(sizeof(Index)*8-1))-1)<<1)+1);
  };

  template <class Cf, class Index>
    const Index imonomial<Cf,Index>::private_max_index_;

/// Indexed polynomial.
  template <class Cf, class Index, class Expo>
    class ipoly
  {
      BOOST_STATIC_ASSERT(boost::is_integral<Index>::value);
      BOOST_STATIC_ASSERT(boost::is_pod<Index>::value);
      BOOST_STATIC_ASSERT(!(boost::integer_traits<Index>::is_signed));
      BOOST_STATIC_ASSERT(boost::is_integral<Expo>::value);
      BOOST_STATIC_ASSERT(boost::is_pod<Expo>::value);
      BOOST_STATIC_ASSERT(!(boost::integer_traits<Expo>::is_signed));
    public:
      typedef std::vector<Expo> vector_expo;
      typedef imonomial<Cf,Index> im_type;
      typedef std::vector<im_type> vector_imonomial;
      ipoly():private_width_(0),private_degree_(0),private_vi_()
        {}
      ipoly(const Cf &value):private_width_(0),private_degree_(0),private_vi_((size_t)1)
      {
        private_vi_[0]=im_type(value,0);
      }
      ipoly(const ipoly &p):private_width_(p.private_width_),private_degree_(p.private_degree_),
        private_vi_(p.private_vi_)
        {}
      ~ipoly()
        {}
// Data members.
    private:
      unsigned short int  private_width_;
      Expo                private_degree_;
      vector_imonomial    private_vi_;
  };
}

#endif
