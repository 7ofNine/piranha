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

#ifndef PIRANHA_TRIG_SPARSE_ARRAY_H
#define PIRANHA_TRIG_SPARSE_ARRAY_H

#include <boost/functional/hash/hash.hpp>
#include <vector>

#include "../../common_typedefs.h"
#include "../../p_assert.h"
#include "../../psymbol.h"
#include "../../utils.h" // lexical_converter.

namespace piranha
{
/// Sparse array trigonometric class.
  class trig_sparse_array
  {
// TODO: try to replace this with struct containing int16[2], in order to employ packing techniques for
// hashing and equality testing.
      typedef std::pair<trig_size_t,int16> pair;
      typedef std::vector<pair> container_type;
      typedef container_type::iterator iterator;
      typedef container_type::const_iterator const_iterator;
    public:
// Start INTERFACE definition.
//-------------------------------------------------------
      typedef int16 value_type;
// Ctors.
/// Default ctor.
      trig_sparse_array():private_flavour_(true) {}
/// Copy ctor.
      trig_sparse_array(const trig_sparse_array &ts):
        private_flavour_(ts.g_flavour()),private_container_(ts.private_container_) {}
      trig_sparse_array(const deque_string &);
      ~trig_sparse_array() {}
// Getters.
      bool &s_flavour()
      {
        return private_flavour_;
      }
      const bool &g_flavour() const
      {
        return private_flavour_;
      }
      int16 at(trig_size_t) const;
      size_t actual_width() const;
// I/O.
      void print_plain(std::ostream &, const vector_psym_p &) const;
      void print_latex(std::ostream &, const vector_psym_p &) const;
// Manip.
/// Assign vector of multipliers.
/**
 * Argument type T must support random access through operator[].
 */
      template <class T>
        void assign_int_vector(const T &v)
      {
        private_container_.resize(0);
        p_assert(empty());
        const size_t w=v.size();
        for (size_t i=0;i<w;++i)
        {
          if (v[i] != 0)
          {
            private_container_.push_back(pair(i,v[i]));
          }
        }
      }
      void append_args(const size_t &) {}
      void prepend_args(const size_t &n)
      {
        const iterator it_f=end();
        for (iterator it=begin();it!=it_f;++it)
        {
// Prepending arguments implies an increase in the index number of each element.
          it->first+=(trig_size_t)n;
        }
      }
      void increase_size(const size_t &) {}
      void invert_sign();
// Probing.
      template <class DerivedPs>
        double density(const DerivedPs &p) const
      {
        switch (p.trig_width())
        {
          case 0:
            return 0;
          default:
            return ((double)private_container_.size())/p.trig_width();
        }
      }
      double freq(const vector_psym_p &) const;
      double phase(const vector_psym_p &) const;
      template <class Series>
        double t_eval(const double &, const Series &) const;
      template <class TrigEvaluator>
        double t_eval(TrigEvaluator &) const;
      short int sign() const;
      size_t hasher() const;
      bool is_zero() const;
      template <class Series>
        bool is_ignorable(const Series &) const;
      bool smaller(const size_t &) const;
      bool larger(const size_t &) const;
      bool is_compatible(const size_t &) const;
      size_t data_footprint() const;
      template <class Series>
        bool checkup(const Series &) const;
      bool operator==(const trig_sparse_array &) const;
      bool operator<(const trig_sparse_array &) const;
// Math.
      void trigmult(const trig_sparse_array &, trig_sparse_array &, trig_sparse_array &) const;
      trig_sparse_array &operator=(const trig_sparse_array &ts2)
      {
        if (&ts2!=this)
        {
          s_flavour()=ts2.g_flavour();
          private_container_=ts2.private_container_;
        }
        return *this;
      }
      trig_sparse_array &operator*=(const int &);
// End INTERFACE definition.
//-------------------------------------------------------
      trig_sparse_array &operator+=(const trig_sparse_array &);
      trig_sparse_array &operator-=(const trig_sparse_array &);
      void dump() const
      {
        const const_iterator it_f=end();
        for (const_iterator it=begin();it!=it_f;++it)
        {
          std::cout << it->first << "," << it->second << '\t';
        }
        std::cout << '\n';
      }
    private:
      const_iterator begin() const
      {
        return private_container_.begin();
      }
      const_iterator end() const
      {
        return private_container_.end();
      }
      iterator begin()
      {
        return private_container_.begin();
      }
      iterator end()
      {
        return private_container_.end();
      }
      bool empty() const
      {
        return private_container_.empty();
      }
      const_iterator find(trig_size_t) const;
      template <class Modifier>
        trig_sparse_array &algebraic_sum(const trig_sparse_array &t2);
// Functors used in generic algebraic addition routine.
      struct sign_modifier_plus
      {
        static const int16 &mod(const int16 &value)
        {
          return value;
        }
        static void mod(const container_type &)
        {}
        static void mod(container_type &c, iterator it1, const_iterator it2, const_iterator it2_f)
        {
          c.insert(it1,it2,it2_f);
        }
      };
      struct sign_modifier_minus
      {
        static int16 mod(const int16 &value)
        {
          return (-value);
        }
        static void mod(container_type &c)
        {
          const size_t w=c.size();
          for (size_t i=0;i<w;++i)
          {
            c[i].second=-c[i].second;
          }
        }
        static void mod(container_type &c, iterator &, const_iterator &it2, const const_iterator &it2_f)
        {
          for (;it2!=it2_f;++it2)
          {
            c.push_back(pair(it2->first,-it2->second));
          }
        }
      };
    private:
      bool            private_flavour_;
      container_type  private_container_;
  };

/// Constructor from piranha::deque_string.
  inline trig_sparse_array::trig_sparse_array(const deque_string &sd):private_flavour_(true),private_container_()
  {
    const size_t w=sd.size();
    if (w==0)
    {
      std::cout << "Warning: constructing empty trig_sparse_array." << std::endl;
      std::abort();
      return;
    }
// TODO: check this.
//    private_container_.reserve(w);
// Now we know  w >= 1.
    int16 tmp_mult;
    for (size_t i=0;i<w-1;++i)
    {
      tmp_mult=utils::lexical_converter<int16>(sd[i]);
      if (tmp_mult != 0)
      {
        private_container_.push_back(pair(i,tmp_mult));
      }
    }
// Take care of flavour.
    if (*sd.back().c_str()=='s')
    {
      s_flavour()=false;
    }
  }

/// Find iterator corresponding to index n. O(n) complexity operation.
  inline trig_sparse_array::const_iterator trig_sparse_array::find(trig_size_t n) const
  {
    const const_iterator it_f=end();
    const_iterator retval=it_f;
    for (const_iterator it=begin();it!=end();++it)
    {
      if (it->first >= n)
      {
// Maybe we found the element. If so assign retval, otherwise...
        if (it->first == n)
        {
          retval=it;
        }
// ... just break, retval was previously assigned to end().
        break;
      }
    }
    return retval;
  }

// Getters implementations.
  inline int16 trig_sparse_array::at(trig_size_t n) const
  {
    const const_iterator it=find(n);
    if (it == end())
    {
      return 0;
    }
    return it->second;
  }

/// Actual width.
  inline size_t trig_sparse_array::actual_width() const
  {
    size_t retval;
    if (private_container_.size() == 0)
    {
      retval=0;
    }
    else
    {
      retval=private_container_.back().first;
    }
    return retval;
  }

// I/O implementations.
// TODO: optimize this, we are calling multiplier which is O(n) but this can be done more efficiently.
  inline void trig_sparse_array::print_plain(std::ostream &out_stream, const vector_psym_p &tv) const
  {
    stream_manager::setup_print(out_stream);
    const_iterator it=begin();
    for (trig_size_t i=0;i<tv.size();++i)
    {
      out_stream << at(i) << stream_manager::data_separator();
    }
    switch (g_flavour())
    {
      case true:
        out_stream << "c";
        break;
      case false:
        out_stream << "s";
    }
  }

  inline void trig_sparse_array::print_latex(std::ostream &out_stream, const vector_psym_p &tv) const
  {
    stream_manager::setup_print(out_stream);
    switch (g_flavour())
    {
      case true:
        out_stream << "c&";
        break;
      case false:
        out_stream << "s&";
    }
    std::string tmp("$");
    const const_iterator it_f=end(), it_b=begin();
    for (const_iterator it=it_b;it!=it_f;++it)
    {
      if (it->second > 0 and it != it_b)
      {
        tmp.append("+");
      }
      if (it->second == -1)
      {
        tmp.append("-");
      }
      else if (it->second == 1) {}
      else
      {
        tmp.append(boost::lexical_cast<std::string>(it->second));
      }
      p_assert(tv.size()>it->first);
      tmp.append(tv[it->first]->name());
    }
    tmp.append("$");
// If we did not write anything erase math markers.
    if (tmp == "$$")
    {
      tmp.clear();
    }
    out_stream << tmp;
  }

  inline void trig_sparse_array::invert_sign()
  {
    const iterator it_f=end();
    for (iterator it=begin();it!=it_f;++it)
    {
      it->second=-it->second;
    }
  }

// Probing implementations.
  inline double trig_sparse_array::freq(const vector_psym_p &v) const
  {
    double retval=0.;
    const const_iterator it_f=end();
    for (const_iterator it=begin();it!=it_f;++it)
    {
// We must be sure that there actually is a freq in every symbol we are going to use.
      if (v[it->first]->poly_eval().size() > 1)
      {
        retval+=it->second*v[it->first]->poly_eval()[1];
      }
    }
    return retval;
  }

  inline double trig_sparse_array::phase(const vector_psym_p &v) const
  {
    double retval=0.;
    const const_iterator it_f=end();
    for (const_iterator it=begin();it!=it_f;++it)
    {
// We must be sure that there actually is a phase in every symbol we are going to use.
      if (v[it->first]->poly_eval().size() > 0)
      {
        retval+=it->second*v[it->first]->poly_eval()[0];
      }
    }
    return retval;
  }

  template <class Series>
    inline double trig_sparse_array::t_eval(const double &t, const Series &s) const
  {
    double retval=0.;
    const const_iterator it_f=end();
    for (const_iterator it=begin();it!=it_f;++it)
    {
      retval+=it->second*s.trig_args()[it->first]->t_eval(t);
    }
    switch (g_flavour())
    {
      case true:
        return std::cos(retval);
      default:
        return std::sin(retval);
    }
  }

  template <class TrigEvaluator>
    inline double trig_sparse_array::t_eval(TrigEvaluator &te) const
  {
    complex_double retval(1.);
    const const_iterator it_f=end();
    for (const_iterator it=begin();it!=it_f;++it)
    {
      retval*=te.request_complexp(it->first,it->second);
    }
    switch (g_flavour())
    {
      case true:
        return retval.real();
      default:
        return retval.imag();
    }
  }

  inline short int trig_sparse_array::sign() const
  {
    short int retval=1;
    if (!empty() && begin()->second < 0)
    {
      retval=-1;
    }
    return retval;
 }

  inline bool trig_sparse_array::operator<(const trig_sparse_array &l2) const
  {
    if (g_flavour() < l2.g_flavour())
    {
      return true;
    }
    else if (g_flavour() > l2.g_flavour())
    {
      return false;
    }
    const const_iterator it_f1=end(), it_f2=l2.end();
    const_iterator it1=begin(), it2=l2.begin();
    while (it1!=it_f1 && it2!=it_f2)
    {
      if (it1->first < it2->first)
      {
        return (it1->second < 0);
      }
      if (it2->first < it1->first)
      {
        return (it2->second > 0);
      }
      if (it1->second < it2->second)
      {
        return true;
      }
      if (it2->second < it1->second)
      {
        return false;
      }
      ++it1;
      ++it2;
    }
    if (it1!=it_f1)
    {
      return (it1->second < 0);
    }
    if (it2!=it_f2)
    {
      return (it2->second > 0);
    }
    return false;
  }

  inline size_t trig_sparse_array::hasher() const
  {
    size_t seed=g_flavour();
    const const_iterator it_f=end();
    for (const_iterator it=begin();it!=it_f;++it)
    {
      boost::hash_combine(seed,it->first);
      boost::hash_combine(seed,it->second);
    }
    return seed;
  }

  inline bool trig_sparse_array::is_zero() const
  {
    return empty();
  }

  template <class Series>
    inline bool trig_sparse_array::is_ignorable(const Series &) const
  {
    return (is_zero() and !g_flavour());
  }

  inline bool trig_sparse_array::smaller(const size_t &) const
  {
    return false;
  }

  inline bool trig_sparse_array::larger(const size_t &) const
  {
    return false;
  }

  inline bool trig_sparse_array::is_compatible(const size_t &) const
  {
    return true;
  }

  inline size_t trig_sparse_array::data_footprint() const
  {
    return (sizeof(pair)*private_container_.size());
  }

  template <class Series>
    inline bool trig_sparse_array::checkup(const Series &s) const
  {
    const const_iterator it_f=end();
// Let's check there are not zero elements.
    for (const_iterator it=begin();it!=it_f;++it)
    {
      if (it->second == 0)
      {
        std::cout << "Zero element in trig_sparse_array." << std::endl;
        return false;
      }
    }
// Now check that we do not have elements in excess.
    if (actual_width() > s.trig_width())
    {
      std::cout << "Too many elements in trig_sparse_array for given size." << std::endl;
      return false;
    }
    return true;
  }

  inline bool trig_sparse_array::operator==(const trig_sparse_array &l2) const
  {
    if (g_flavour() != l2.g_flavour())
    {
      return false;
    }
    return (private_container_ == l2.private_container_);
  }

  inline trig_sparse_array &trig_sparse_array::operator*=(const int &n)
  {
    if (n==0)
    {
      private_container_.clear();
    }
    else
    {
      const iterator it_f=end();
      for (iterator it=begin();it!=it_f;++it)
      {
        it->second*=n;
      }
    }
    return *this;
  }

  template <class Modifier>
    inline trig_sparse_array &trig_sparse_array::algebraic_sum(const trig_sparse_array &t2)
  {
    p_assert(this!=&t2);
    if (t2.empty())
    {
      return *this;
    }
    if (empty())
    {
      private_container_=t2.private_container_;
      Modifier::mod(private_container_);
      return *this;
    }
    const size_t w1=private_container_.size(), w2=t2.private_container_.size();
// At this point both vectors have non-zero size. Reserve in advance, in most cases we will be adding elements.
    private_container_.reserve(std::max(w1,w2));
    iterator it1=begin();
    const_iterator it2=t2.begin();
    const const_iterator it2_f=t2.end();
    while (it2 != it2_f)
    {
// We are at the end of this, insert the remaining elements and bail out of the cycle.
      if (it1 == end())
      {
        Modifier::mod(private_container_,it1,it2,it2_f);
        break;
      }
      else
      {
// Same index, add/subtract.
        if (it1->first == it2->first)
        {
          it1->second+=Modifier::mod(it2->second);
// If we modified to zero, we have to destroy the element.
          if (it1->second == 0)
          {
            it1=private_container_.erase(it1); // it1 now points to the element after the erased one (the latter half of the
                                              // vector was moved down by one position).
          }
          else
// We performed a non-destructive modification. Increase it1.
          {
            ++it1;
          }
// Increase it2, it was added to this.
          ++it2;
        }
// There is an element which is not present in this and which goes before it1.
        else if (it1->first > it2->first)
        {
// it1 will point to the newly inserted element.
          it1=private_container_.insert(it1,pair(it2->first,Modifier::mod(it2->second)));
          ++it1;
          ++it2;
        }
// it2's index is after it1's. We don't do anything since we don't know if next elements of t2 will be packed or inserted.
        else /* if (it1->first < it2->first) */
        {
          ++it1;
        }
      }
    }
    return *this;
/*    while (it2!=it_f2)
    {
// We are at the end of this, insert new term and do not increase iterators on this.
      if (it1==end())
      {
        old_it1=insert_after(it2->first,Modifier::mod(it2->second),old_it1);
        ++it2;
      }
// We found a duplicate, modify it.
      else if (it1->first==it2->first)
      {
        old_it1=modify(it1->second+Modifier::mod(it2->second),old_it1,it1);
// Take care in case we erased the element when modifying, we may have run out of elements.
        if (old_it1==end())
        {
// We erased the first element of this' list.
          it1=begin();
        }
        else
        {
// We modified to a non-zero value an element of this' list or erased an element else than the first one.
          it1=old_it1;
          ++it1;
        }
        ++it2;
      }
// this has gone past t2.
      else if (it1->first > it2->first)
      {
        old_it1=insert_after(it2->first,Modifier::mod(it2->second),old_it1);
        ++it2;
      }
// t2 has gone past this.
      else if (it1->first < it2->first)
      {
        old_it1=it1;
        ++it1;
      }
    }
    return *this;*/
  }

  inline trig_sparse_array &trig_sparse_array::operator+=(const trig_sparse_array &t2)
  {
    return algebraic_sum<sign_modifier_plus>(t2);
  }

  inline trig_sparse_array &trig_sparse_array::operator-=(const trig_sparse_array &t2)
  {
    return algebraic_sum<sign_modifier_minus>(t2);
  }

  inline void trig_sparse_array::trigmult(const trig_sparse_array &l2, trig_sparse_array &ret1, trig_sparse_array &ret2) const
  {
    ret1=*this;
    ret1-=l2;
    ret2=*this;
    ret2+=l2;
  }

// NOTE: maybe we can place this struct inside trig_sparse_array and require that all trig_args provide
// it, as an interface request? But then all classes that are expected to use the hasher should then
// multi-inherit from a base hasher class? Otherwise how can we overload hash-value? Mh, maybe we
// can specificy the hasher in multiindex containers. Investigate.
/// Overload of hash_value function for piranha::trig_sparse_array.
/**
 * To be used in piranha::base_pseries for the hashed index.
 */
  inline size_t hash_value(const trig_sparse_array &tm)
  {
    return tm.hasher();
  }
}
#endif                                            // PIRANHA_trig_sparse_array_H
