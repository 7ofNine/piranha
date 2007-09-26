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

#include <algorithm>
#include <boost/integer_traits.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/identity.hpp>
#include <boost/multi_index_container.hpp>
#include <boost/static_assert.hpp>
#include <boost/type_traits/is_integral.hpp>
#include <boost/type_traits/is_pod.hpp>
#include <cmath>
#include <ext/pool_allocator.h>
#include <limits>
#include <vector>

#include "p_assert.h"

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
      Cf &s_cf()
      {
        return private_cf_;
      }
      static const Index &g_max_index()
      {
        return private_max_index_;
      }
    private:
      imonomial()
        {}
    private:
      Index               private_index_;
      Cf                  private_cf_;
      static const Index  private_max_index_ = ((((((Index)1)<<(sizeof(Index)*8-1))-1)<<1)+1);
  };

  template <class Cf, class Index>
    const Index imonomial<Cf,Index>::private_max_index_;

/// Lightweight monomial class used during polynomial multiplication.
  template <class Cf, class Index>
    struct mutable_im
  {
      mutable_im(const Index &i, const Cf &c):first(i),second(c)
        {}
      bool operator==(const mutable_im &m) const
      {
        return (first == m.first);
      }
      bool operator<(const mutable_im &m) const
      {
        return (first < m.first);
      }
      Index       first;
      mutable Cf  second;
    private:
// Make it private so we know it is never called.
      mutable_im()
        {}
  };

  template <class Cf, class Index>
    inline const Index &hash_value(const mutable_im<Cf,Index> &m)
  {
    return m.first;
  }

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
        private_vi_()
      {
        p_assert(v.size() <= USHRT_MAX);
        const usint w=g_width();
        for (usint i=0;i<w;++i)
        {
          private_degree_+=v[i];
        }
//std::cout << "Assigned with properties: " << private_width_ << '\t' << private_degree_ << '\n';
        Index tmp_index;
        encode(v,tmp_index);
        private_vi_.push_back(im_type(value,tmp_index));
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
          private_width_=p.g_width();
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
      bool empty() const
      {
        return private_vi_.empty();
      }
      void clear()
      {
        private_vi_.clear();
        private_degree_=0;
      }
      size_t g_length() const
      {
        return private_vi_.size();
      }
      const usint &g_width() const
      {
        return private_width_;
      }
      const usint &g_degree() const
      {
        return private_degree_;
      }
      void decode(const Index &code, vector_expo &v) const
      {
        p_assert(v.size() == g_width());
        if (private_width_ == 0)
        {
// No need to do anything, v is zero-dimensioned.
          return;
        }
        const Expo max_d=private_max_n_cache_.g_max_n(private_width_);
//std::cout << "decoding " << code << '\n';
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
      ipoly &operator+=(const ipoly &p)
      {
        algebraic_sum<sign_modifier_plus>(p);
        return *this;
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
      void refresh_degree()
      {
        vector_expo v((size_t)private_width_);
        const iterator it_f=end();
        private_degree_=0;
        Expo candidate;
        usint i;
        for (iterator it=begin();it!=it_f;++it)
        {
          decode(it->g_index(),v);
          candidate=0;
          for (i=0;i<private_width_;++i)
          {
            candidate+=v[i];
          }
          if (candidate > private_degree_)
          {
            private_degree_=candidate;
          }
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
//std::cout << "Max representable degree is: " << max_d << '\n';
        retval=0;
        for (usint i=0;i<private_width_;++i)
        {
          retval+=v[i]*integral_npow_cache<Index>::request(i,max_d+1);
        }
//std::cout << "encoded to " << retval << '\n';
      }
// Boilerplate for algebraic sum.
      struct sign_modifier_plus
      {
        static const Cf &mod(const Cf &value)
        {
          return value;
        }
        static void mod(const vector_imonomial &)
        {}
        static void mod(vector_imonomial &v, iterator it1, const_iterator it2, const_iterator it2_f)
        {
          v.insert(it1,it2,it2_f);
        }
      };
      struct sign_modifier_minus
      {
        static Cf mod(const Cf &value)
        {
          return (-value);
        }
        static void mod(vector_imonomial &v)
        {
          const iterator it_f=v.end();
          for (iterator it=v.begin();it!=it_f;++it)
          {
            it->second=-(it->second);
          }
        }
        static void mod(vector_imonomial &v, iterator &, const_iterator &it2, const const_iterator &it2_f)
        {
          for (;it2!=it2_f;++it2)
          {
            v.push_back(im_type(it2->g_cf(),-it2->g_index()));
          }
        }
      };
      template <class Modifier>
        void algebraic_sum(const ipoly &p)
      {
        if (p.empty())
        {
          return;
        }
        if (empty())
        {
          private_vi_=p.private_vi_;
          Modifier::mod(private_vi_);
          return;
        }
        if (this == &p)
        {
          ipoly tmp(p);
          algebraic_sum<Modifier>(tmp);
          return;
        }
        const size_t w1=private_vi_.size(), w2=p.private_vi_.size();
// TODO: check about this, maybe we risk allocating too much after repeated additions?
// At this point both vectors have non-zero size. Reserve in advance, in most cases we will be adding elements.
        private_vi_.reserve(std::max(w1,w2));
        iterator it1=begin();
        const_iterator it2=p.begin();
        const const_iterator it2_f=p.end();
        while (it2 != it2_f)
        {
// We are at the end of this, insert the remaining elements and bail out of the cycle.
          if (it1 == end())
          {
            Modifier::mod(private_vi_,it1,it2,it2_f);
            break;
          }
          else
          {
// Same key, add/subtract.
            if (it1->g_index() == it2->g_index())
            {
              it1->s_cf()+=Modifier::mod(it2->g_cf());
// If we modified to zero, we have to destroy the element.
// FIXME: replace with numerical zero.
              if (it1->g_cf() == 0)
              {
                it1=private_vi_.erase(it1); // it1 now points to the element after the erased one (the latter half of the
                                            // vector was moved down by one position).
              }
              else
// We performed a non-destructive modification. Increase it1.
              {
                ++it1;
              }
// Increase ip, it was added to this.
              ++it2;
            }
// There is an element which is not present in this and which goes before it1.
            else if (it1->g_index() > it2->g_index())
            {
// it1 will point to the newly inserted element.
              it1=private_vi_.insert(it1,im_type(Modifier::mod(it2->g_cf()),it2->g_index()));
              ++it1;
              ++it2;
            }
// ip's index is after it1's. We don't do anything since we don't know if next elements of p will be packed or inserted.
            else /* if (it1->first < ip->first) */
            {
              ++it1;
            }
          }
        }
// Refresh degree: we may have introduced higher degree monomials, and we may have destroyed others.
// Just recalc it explicitly from scracth.
        refresh_degree();
      }
// Multiplication boilerplate.
      struct index_sorter
      {
        bool operator()(const im_type &m1, const im_type &m2) const
        {
          return (m1.g_index() < m2.g_index());
        }
      };
      typedef boost::multi_index_container<
        mutable_im<Cf,Index>,
        boost::multi_index::indexed_by<
          boost::multi_index::hashed_unique<boost::multi_index::identity<mutable_im<Cf,Index> > >
        >,
      __gnu_cxx::__pool_alloc<mutable_im<Cf,Index> > > hash_map;
      void mult_by(const ipoly &p)
      {
        if (empty())
        {
          p_assert(private_degree_ == 0);
          return;
        }
        if (p.empty())
        {
          p_assert(p.g_degree() == 0);
          clear();
          return;
        }
        if (this == &p)
        {
          ipoly tmp(*this);
          mult_by(tmp);
          return;
        }
        p_assert(private_width_ == p.private_width_);
        const Expo new_degree=g_degree()+p.g_degree();
//        std::cout << "New degree will be: " << new_degree << '\n';
        if (new_degree > private_max_n_cache_.g_max_n(private_width_))
        {
          std::cout << "FATAL: polynomial multiplication results in overflow degree." << std::endl;
          std::abort();
        }
        hash_map tmp;
        typedef typename hash_map::iterator hm_it;
        std::pair<hm_it,bool> insert_res;
        const const_iterator it1_f=end(), it2_f=p.end();
        for (const_iterator it1=begin();it1!=it1_f;++it1)
        {
          for (const_iterator it2=p.begin();it2!=it2_f;++it2)
          {
            mutable_im<Cf,Index> tmp_m(it1->g_index()+it2->g_index(),it1->g_cf()*it2->g_cf());
            insert_res=tmp.insert(tmp_m);
            if (!(insert_res.second))
            {
// We found a duplicate.
              insert_res.first->second+=tmp_m.second;
            }
          }
        }
// Now place the result in a return value.
        ipoly retval;
        retval.private_width_=private_width_;
        retval.private_vi_.reserve(tmp.size());
        const hm_it it_f=tmp.end();
        for (hm_it it=tmp.begin();it!=it_f;++it)
        {
          retval.private_vi_.push_back(im_type(it->second,it->first));
        }
        retval.private_degree_=new_degree;
// Sort result according to index.
        std::sort(retval.begin(),retval.end(),index_sorter());
        swap(retval);
std::cout << "Length: " << g_length() << '\n';
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
//std::cout << "Max expo is: " << private_max_expo_ << '\n';
//std::cout << private_container_[0] << ',';
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
//std::cout << private_container_[i] << ',';
          }
//std::cout << '\n';
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
      /*const*/ usint               private_width_;
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
    typedef typename piranha::ipoly<Cf,Index,Expo>::const_iterator const_iterator;
    typename piranha::ipoly<Cf,Index,Expo>::vector_expo ev(p.g_width());
    const const_iterator it_f=p.end();
    for (const_iterator it=p.begin();it!=it_f;++it)
    {
      std::cout << it->g_cf() << "\t|\t";
      p.decode(it->g_index(),ev);
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
