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
#include <cmath>
#include <limits>
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
        p_assert(arg >= 0);
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
// Getters.
      const Cf &g_cf() const
      {
        return private_cf_;
      }
      const Index &g_index() const
      {
        return private_index_;
      }
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
      typedef unsigned short int usint;
      typedef typename vector_imonomial::iterator iterator;
      typedef typename vector_imonomial::const_iterator const_iterator;
      ipoly():private_width_(0),private_degree_(0),private_vi_()
        {}
      ipoly(const Cf &value):private_width_(0),private_degree_(0),private_vi_((size_t)1)
      {
        private_vi_[0]=im_type(value,0);
      }
      ipoly(const Cf &value, const vector_expo &v):private_width_(v.size()),private_degree_(0),
        private_vi_((size_t)1)
      {
        p_assert(v.size() <= USHRT_MAX);
        for (usint i=0;i<private_width_;++i)
        {
          private_degree_+=v[i];
        }
std::cout << "Assigned with properties: " << private_width_ << '\t' << private_degree_ << '\n';
        Index tmp_index;
        encode(v,tmp_index);
        im_type tmp_im(value,tmp_index);
        private_vi_[0]=tmp_im;
        degree_check();
      }
      ipoly(const ipoly &p):private_width_(p.private_width_),private_degree_(p.private_degree_),
        private_vi_(p.private_vi_)
        {}
      ~ipoly()
        {}
      ipoly &operator=(const ipoly &p)
      {
        if (this != &p)
        {
          private_width_=p.private_width_;
          private_degree_=p.private_degree_;
          private_vi_=p.private_vi_;
        }
        return *this;
      }
      iterator begin()
      {
        return private_vi_.begin();
      }
      iterator end()
      {
        return private_vi_.end();
      }
      const_iterator begin() const
      {
        return private_vi_.begin();
      }
      const_iterator end() const
      {
        return private_vi_.end();
      }
      const usint &g_width() const
      {
        return private_width_;
      }
      size_t g_length() const
      {
        return private_vi_.size();
      }
      const im_type &operator[](const size_t &n) const
      {
        return private_vi_[n];
      }
      void decode(const Index &code, vector_expo &v) const
      {
        p_assert(v.size() == private_width_);
        if (private_width_ == 0)
        {
          return;
        }
        const Expo max_d=private_max_n_cache_.g_max_n(private_width_);
std::cout << "decoding " << code << '\n';
        for (usint i=0;i<private_width_-1;++i)
        {
          v[i]=(code%integral_npow_cache<Index>::request(i+1,max_d+1))/
            integral_npow_cache<Index>::request(i,max_d+1);
        }
        v[private_width_-1]=code/(integral_npow_cache<Index>::request(private_width_-1,max_d+1));
/*std::cout << "decoded to: ";
for (size_t i=0;i<v.size();++i)
{
  std::cout << v[i] << ',';
}
std::cout << '\n';*/
      }
      ipoly &operator*=(const ipoly &p)
      {
        mult_by(p);
        return *this;
      }
    private:
      void degree_check() const
      {
        if (private_degree_ > private_max_n_cache_.g_max_n(private_width_))
        {
          std::cout << "FATAL: polynomial degree is too large." << std::endl;
          std::abort();
        }
      }
      ipoly &swap(ipoly &p)
      {
        if (this != &p)
        {
          std::swap(private_width_,p.private_width_);
          std::swap(private_degree_,p.private_degree_);
          private_vi_.swap(p.private_vi_);
        }
        return *this;
      }
      void encode(const vector_expo &v, Index &retval) const
      {
        p_assert(v.size() == private_width_);
// Maximum representable degree.
        const Expo max_d=private_max_n_cache_.g_max_n(private_width_);
std::cout << "Max representable degree is: " << max_d << '\n';
        retval=0;
        for (usint i=0;i<private_width_;++i)
        {
          retval+=v[i]*integral_npow_cache<Index>::request(i,max_d+1);
        }
std::cout << "encoded to " << retval << '\n';
      }
      void mult_by(const ipoly &p)
      {
        const const_iterator it1_f=end(), it2_f=p.end();
      }
    class max_n_cache
    {
      public:
        max_n_cache():private_container_(std::floor(std::log(im_type::g_max_index())/std::log(2))+1)
        {
          const size_t w=private_container_.size();
          p_assert(w >= 1);
// With zero-variables polynomials we can go to whatever degree we want.
          private_container_[0]=private_max_expo_;
std::cout << "Max expo is: " << private_max_expo_ << '\n';
std::cout << private_container_[0] << ',';
          double tmp;
          for (size_t i=1;i<w;++i)
          {
            tmp=std::floor(std::pow(im_type::g_max_index(),1./i))-1;
            if (tmp > private_max_expo_)
            {
              private_container_[i]=private_max_expo_;
            }
            else
            {
              private_container_[i]=(Expo)(tmp);
            }
std::cout << private_container_[i] << ',';
          }
std::cout << '\n';
        }
        const Expo &g_max_n(const usint &n) const
        {
          p_assert(n < private_container_.size());
          return private_container_[n];
        }
      private:
        vector_expo private_container_;
    };
// Data members.
    private:
      const usint               private_width_;
      Expo                      private_degree_;
      vector_imonomial          private_vi_;
      static const max_n_cache  private_max_n_cache_;
      static const Expo         private_max_expo_ = ((((((Expo)1)<<(sizeof(Expo)*8-1))-1)<<1)+1);
  };

  template <class Cf, class Index, class Expo>
    const typename ipoly<Cf,Index,Expo>::max_n_cache ipoly<Cf,Index,Expo>::private_max_n_cache_;

  template <class Cf, class Index, class Expo>
    const Expo ipoly<Cf,Index,Expo>::private_max_expo_;
}

namespace std
{
  template <class Cf, class Index, class Expo>
    ostream &operator<<(ostream &s, const piranha::ipoly<Cf,Index,Expo> &p)
  {
    typename piranha::ipoly<Cf,Index,Expo>::vector_expo ev(p.g_width());
    for (size_t i=0;i<p.g_length();++i)
    {
      std::cout << p[i].g_cf() << "\t|\t";
      p.decode(p[i].g_index(),ev);
      for (size_t j=0;j<ev.size();++j)
      {
        std::cout << ev[ev.size()-j-1] << '\t';
      }
      std::cout << std::endl;
    }
    return s;
  }
}

#endif
