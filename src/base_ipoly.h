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

#ifndef PIRANHA_BASE_IPOLY_H
#define PIRANHA_BASE_IPOLY_H

#include <boost/integer_traits.hpp>
#include <boost/static_assert.hpp>
#include <boost/type_traits/is_integral.hpp>
#include <boost/type_traits/is_pod.hpp>
#include <ext/hash_set>
#include <limits>
#include <valarray>

#include "imonomial.h"
#include "integral_npow_cache.h"
#include "p_assert.h"
#include "utils.h"

namespace piranha
{
/// Base indexed polynomial.
  template <class Cf, class Index, class Expo, class Derived>
    class base_ipoly
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
      typedef __gnu_cxx::hash_set<im_type,monomial_hasher<im_type> > container_type;
      typedef unsigned short int usint;
      typedef typename container_type::iterator iterator;
      typedef typename container_type::const_iterator const_iterator;
      base_ipoly():private_degree_(0),private_container_()
        {}
      base_ipoly(const Cf &value):private_degree_(0),private_container_()
      {
        private_container_.insert(im_type(value,0));
      }
      base_ipoly(const base_ipoly &p):private_degree_(p.private_degree_),private_container_(p.private_container_)
        {}
      base_ipoly(const Derived &p):private_degree_(p.private_degree_),private_container_(p.private_container_)
        {}
      ~base_ipoly()
        {}
      base_ipoly &operator=(const base_ipoly &p)
      {
        if (this != &p)
        {
          common_assignment(p);
        }
        return *this;
      }
      const_iterator begin() const
      {
        return private_container_.begin();
      }
      const_iterator end() const
      {
        return private_container_.end();
      }
      bool empty() const
      {
        return private_container_.empty();
      }
      void clear()
      {
        private_container_.clear();
        private_degree_=0;
      }
      size_t length() const
      {
        return private_container_.size();
      }
      const usint &g_degree() const
      {
        return private_degree_;
      }
      void decode(const Index &code, vector_expo &v) const
      {
        const usint w=static_cast<Derived const *>(this)->g_width();
        p_assert(v.size() == w);
        if (w == 0)
        {
// No need to do anything, v is zero-dimensioned.
          return;
        }
        const Expo max_d=private_max_n_cache_.g_max_n(w);
//std::cout << "decoding " << code << '\n';
        for (usint i=0;i<w-1;++i)
        {
          v[i]=(code%integral_npow_cache<Index>::request(i+1,max_d+1))/
            integral_npow_cache<Index>::request(i,max_d+1);
        }
        v[w-1]=code/(integral_npow_cache<Index>::request(w-1,max_d+1));
/*std::cout << "decoded to: ";
for (size_t i=0;i<v.size();++i)
{
  std::cout << v[i] << ',';
}
std::cout << '\n';*/
      }
    protected:
      void builder_from_vector(const Cf &value, const vector_expo &v)
      {
        if (v.size() > USHRT_MAX)
        {
          std::cout << "Fatal: cannot build polynomial, vector is too large." << std::endl;
          std::cout << "Please use a vector whose size is smaller than " << USHRT_MAX << std::endl;
          std::abort();
        }
        private_degree_=0;
        const usint w=static_cast<Derived const *>(this)->g_width();
        for (usint i=0;i<w;++i)
        {
          private_degree_+=v[i];
        }
//std::cout << "Assigned with properties: " << private_width_ << '\t' << private_degree_ << '\n';
        Index tmp_index;
        encode(v,tmp_index);
        private_container_.insert(im_type(value,tmp_index));
        degree_check();
      }
      iterator begin()
      {
        return private_container_.begin();
      }
      iterator end()
      {
        return private_container_.end();
      }
      void common_assignment(const Derived &p)
      {
          private_degree_=p.private_degree_;
          private_container_=p.private_container_;
      }
      void degree_check() const
      {
        if (private_degree_ > private_max_n_cache_.g_max_n(static_cast<Derived const *>(this)->g_width()))
        {
          std::cout << "FATAL: polynomial degree is too large." << std::endl;
          std::abort();
        }
      }
      void refresh_degree()
      {
        const usint w=static_cast<Derived const *>(this)->g_width();
        vector_expo v((size_t)w);
        const iterator it_f=end();
        private_degree_=0;
        Expo candidate;
        usint i;
        for (iterator it=begin();it!=it_f;++it)
        {
          decode(it->g_index(),v);
          candidate=0;
          for (i=0;i<w;++i)
          {
            candidate+=v[i];
          }
          if (candidate > private_degree_)
          {
            private_degree_=candidate;
          }
        }
      }
      template <class Modifier>
        void insert(const im_type &m)
      {
// TODO: think about zero detection here.
        if (m.cf == 0)
        {
          return;
        }
        iterator it=private_container_.find(m);
        if (it == end())
        {
// Not a duplicate, insert.
          action_assert(private_container_.insert(im_type(Modifier::mod(m.cf),m.index)).second);
        }
        else
        {
// Duplicate: merge with existing element.
          it->cf+=Modifier::mod(m.cf);
// If the result is zero, erase.
          if (it->cf == 0)
          {
            private_container_.erase(it);
          }
        }
      }
      void encode(const vector_expo &v, Index &retval) const
      {
        const usint w=static_cast<Derived const *>(this)->g_width();
        p_assert(v.size() == w);
// Maximum representable degree.
        const Expo max_d=private_max_n_cache_.g_max_n(w);
//std::cout << "Max representable degree is: " << max_d << '\n';
        retval=0;
        for (usint i=0;i<w;++i)
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
        static void mod(const container_type &)
        {}
      };
      struct sign_modifier_minus
      {
        static Cf mod(const Cf &value)
        {
          return (-value);
        }
        static void mod(container_type &c)
        {
          const iterator it_f=c.end();
          for (iterator it=c.begin();it!=it_f;++it)
          {
            it->cf=-(it->cf);
          }
        }
      };
      template <class Modifier>
        void algebraic_sum(const Derived &p)
      {
        if (p.empty())
        {
          return;
        }
        if (empty())
        {
          private_container_=p.private_container_;
          Modifier::mod(private_container_);
          return;
        }
        if (this == &p)
        {
          Derived tmp(p);
          algebraic_sum<Modifier>(tmp);
          return;
        }
        const const_iterator it2_f=p.end();
        for (const_iterator it2=p.begin();it2!=it2_f;++it2)
        {
          insert<Modifier>(*it2);
        }
// Refresh degree: we may have introduced higher degree monomials, and we may have destroyed others.
// Just recalc it explicitly from scracth.
        refresh_degree();
      }
      void addition(const Derived &p)
      {
        algebraic_sum<sign_modifier_plus>(p);
      }
      void subtraction(const Derived &p)
      {
        algebraic_sum<sign_modifier_minus>(p);
      }
      void mult_by(const Derived &p)
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
          Derived tmp(*static_cast<const Derived *>(this));
          mult_by(tmp);
          return;
        }
        const usint w=static_cast<Derived const *>(this)->g_width();
        p_assert(w == p.g_width());
//        std::cout << "New degree will be: " << new_degree << '\n';
        if (g_degree()+p.g_degree() > private_max_n_cache_.g_max_n(w))
        {
          std::cout << "FATAL: polynomial multiplication results in overflow degree." << std::endl;
          std::abort();
        }
// This ipoly acts just as a container for the private_container_, which will be swapped in at the end of the cycle.
        Derived tmp;
        std::valarray<const_iterator> v_it2;
        utils::array_iter(p,v_it2);
        size_t i;
        const size_t l2=p.length();
        const const_iterator it1_f=end();
        for (const_iterator it1=begin();it1!=it1_f;++it1)
        {
          for (i=0;i<l2;++i)
          {
            tmp.insert<sign_modifier_plus>(im_type(it1->g_cf()*v_it2[i]->g_cf(),it1->g_index()+v_it2[i]->g_index()));
          }
        }
        private_container_.swap(tmp.private_container_);
// Refresh degree: we may have introduced higher degree monomials, and we may have destroyed others.
// Just recalc it explicitly from scracth.
        refresh_degree();
std::cout << "Length: " << length() << '\n';
      }
    class max_n_cache
    {
      public:
        max_n_cache():private_container_(std::floor(std::log(im_type::max_index)/std::log(2))+1)
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
            tmp=std::floor(std::pow(im_type::max_index,1./i))-1;
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
    protected:
      Expo                      private_degree_;
      container_type            private_container_;
      static const max_n_cache  private_max_n_cache_;
      static const Expo         private_max_expo_ = ((((((Expo)1)<<(sizeof(Expo)*8-1))-1)<<1)+1);
  };

  template <class Cf, class Index, class Expo, class Derived>
    const typename base_ipoly<Cf,Index,Expo,Derived>::max_n_cache
    base_ipoly<Cf,Index,Expo,Derived>::private_max_n_cache_;

  template <class Cf, class Index, class Expo, class Derived>
    const Expo base_ipoly<Cf,Index,Expo,Derived>::private_max_expo_;
}

namespace std
{
  template <class Cf, class Index, class Expo, class Derived>
    ostream &operator<<(ostream &s, const piranha::base_ipoly<Cf,Index,Expo,Derived> &p)
  {
    typedef typename piranha::base_ipoly<Cf,Index,Expo,Derived>::const_iterator const_iterator;
    typename piranha::base_ipoly<Cf,Index,Expo,Derived>::vector_expo
      ev(static_cast<Derived const *>(&p)->g_width());
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
