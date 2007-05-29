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

#ifndef PIRANHA_TRIG_SLIST_H
#define PIRANHA_TRIG_SLIST_H

#include <boost/functional/hash/hash.hpp>
#include <ext/slist>

#include "common_typedefs.h"
#include "p_assert.h"
#include "psymbol_manager.h"
#include "utils.h"            // lexical_converter.

namespace piranha
  {
  /// Linear combination of trigonometric arguments. Singly linked list version.
  /**
   * Implementation of the trigonometric arguments interface implemented with a singly linked list.
   * The interface is in common with piranha::trig_array.
   * Multipliers of trigonometric arguments are kept ordered according to their index values.
   * @see piranha::trig_array for the documentation of the interface.
   */
  class trig_slist
    {
    public:
      // Start INTERFACE definition.
      //-------------------------------------------------------
      // Ctors.
      trig_slist()
      {}
      trig_slist(const trig_slist &ts):container_(ts.container_)
      {}
      trig_slist(const deque_string &);
      ~trig_slist()
      {}
      // Getters.
      mult_t multiplier(trig_size_t) const;
      size_t actual_width() const;
      // I/O.
      void print_plain(std::ostream &, const vector_psym_p &) const;
      void print_latex(std::ostream &, const vector_psym_p &) const;
      // Manip.
      // FIXME: move to protected once it is deprecated?
      void insert(trig_size_t,mult_t);
      void add_arg()
      {}
      void resize(const size_t &)
      {}
      void invert_sign();
      // Probing.
      double freq(const vector_psym_p &) const;
      double phase(const vector_psym_p &) const;
      double t_eval(const double &, const vector_psym_p &) const;
      short int sign() const;
      size_t hasher() const;
      bool is_zero() const;
      bool smaller(const size_t &) const;
      bool larger(const size_t &) const;
      bool compatible(const size_t &) const;
      size_t data_footprint() const;
      bool checkup(const size_t &) const;
      bool operator==(const trig_slist &) const;
      bool operator>(const trig_slist &) const;
      bool operator<(const trig_slist &) const;
      // Math.
      void trigmult(const trig_slist &, trig_slist &, trig_slist &) const;
      trig_slist &operator=(const trig_slist &ts2)
      {
        if (&ts2!=this)
          {
            container_=ts2.container_;
          }
        return *this;
      }
      trig_slist &operator*=(mult_t);
      // End INTERFACE definition.
      //-------------------------------------------------------
      trig_slist &operator+=(const trig_slist &);
      trig_slist &operator-=(const trig_slist &);
    private:
      typedef std::pair<trig_size_t,mult_t> pair;
      typedef __gnu_cxx::slist<pair>::iterator iterator;
      typedef __gnu_cxx::slist<pair>::const_iterator const_iterator;
      typedef std::pair<iterator,iterator> find_result;
      typedef std::pair<const_iterator,const_iterator> const_find_result;
      const_iterator begin() const
        {
          return container_.begin();
        }
      const_iterator end() const
        {
          return container_.end();
        }
      iterator begin()
      {
        return container_.begin();
      }
      iterator end()
      {
        return container_.end();
      }
      bool empty() const
        {
          return container_.empty();
        }
      void find(trig_size_t, const_find_result &) const;
      void find(trig_size_t, find_result &);
      void insert_new(trig_size_t, mult_t, const find_result &);
      iterator insert_after(trig_size_t, mult_t, iterator);
      iterator modify(mult_t, iterator, iterator);
      template <short int Sign>
      trig_slist &algebraic_sum(const trig_slist &t2);
    private:
      __gnu_cxx::slist<pair>    container_;
    };


  inline trig_slist::trig_slist(const deque_string &sd)
  {
    const size_t w=sd.size();
    for (size_t i=0;i<w;++i)
      {
        insert(i,utils::lexical_converter<mult_t>(sd[i]));
      }
  }


  /// Find iterator corresponding to index n. O(n) time operation.
  inline void trig_slist::find(trig_size_t n, const_find_result &fr) const
    {
      fr.first=end();
      fr.second=begin();
      const const_iterator it_f=end();
      for (;fr.second!=it_f;++fr.second)
        {
          if (fr.second->first>=n)
            {
              // Maybe we found the element. If so return, otherwise...
              if (fr.second->first>n)
                {
                  // ... element was not found, assign end() to second return value.
                  fr.second=it_f;
                }
              break;
            }
          fr.first=fr.second;
        }
    }


  inline void trig_slist::find(trig_size_t n, find_result &fr)
  {
    fr.first=end();
    fr.second=begin();
    const iterator it_f=end();
    for (;fr.second!=it_f;++fr.second)
      {
        if (fr.second->first>=n)
          {
            // Maybe we found the element. If so return, otherwise...
            if (fr.second->first>n)
              {
                // ... element was not found, assign end() to second return value.
                fr.second=it_f;
              }
            break;
          }
        fr.first=fr.second;
      }
  }


  // Getters implementations.
  inline mult_t trig_slist::multiplier(trig_size_t n) const
    {
      const_find_result fr;
      find(n,fr);
      if (fr.second==end())
        {
          return 0;
        }
      return fr.second->second;
    }


  inline size_t trig_slist::actual_width() const
    {
      const const_iterator it_f=end();
      const_iterator it_prev=it_f;
      for (const_iterator it=begin();it!=it_f;++it)
        {
          if (it_prev==it_f)
            {
              it_prev=it;
            }
          else
            {
              ++it_prev;
            }
        }
      if (it_prev==it_f)
        {
          return 0;
        }
      else
        {
          return it_prev->first;
        }
    }


  // I/O implementations.
  inline void trig_slist::print_plain(std::ostream &out_stream, const vector_psym_p &tv) const
    {
      stream_manager::setup_print(out_stream);
      const_iterator it=begin();
      for (trig_size_t i=0;i<tv.size();++i)
        {
          out_stream << multiplier(i) << stream_manager::data_separator();
        }
    }


  inline void trig_slist::print_latex(std::ostream &out_stream, const vector_psym_p &tv) const
    {
      stream_manager::setup_print(out_stream);
      std::string tmp("$");
      const const_iterator it_f=end(), it_b=begin();
      for (const_iterator it=it_b;it!=it_f;++it)
        {
          if (it->second>0 && it!=it_b)
            {
              tmp.append("+");
            }
          if (it->second==-1)
            {
              tmp.append("-");
            }
          else if (it->second==1)
          {}
          else
            {
              tmp.append(boost::lexical_cast<std::string>(it->second));
            }
          // Add a check.
          p_assert(tv.size()>it->first);
          tmp.append(tv[it->first]->name());
        }
      tmp.append("$");
      // If we did not write anything erase math markers.
      if (tmp=="$$")
        {
          tmp.clear();
        }
      out_stream << tmp;
    }


  // Manipulation implementations.
  inline void trig_slist::insert_new(trig_size_t index, mult_t multiplier, const find_result &fr)
  {
    // Do not insert a new zero multiplier.
    if (multiplier==0)
      {
        return;
      }
    if (fr.first==end())
      {
        container_.push_front(pair(index,multiplier));
      }
    else
      {
        container_.insert_after(fr.first,pair(index,multiplier));
      }
  }


  inline trig_slist::iterator trig_slist::insert_after(trig_size_t index, mult_t multiplier, iterator it)
  {
    // Do not insert a new zero multiplier.
    if (multiplier==0)
      {
        return it;
      }
    if (it==end())
      {
        container_.push_front(pair(index,multiplier));
        return container_.begin();
      }
    else
      {
        return container_.insert_after(it,pair(index,multiplier));
      }
  }


  inline trig_slist::iterator trig_slist::modify(mult_t multiplier, iterator prev, iterator cur)
  {
    if (multiplier!=0)
      {
        cur->second=multiplier;
        return cur;
      }
    else
      {
        if (prev==end())
          {
            p_assert(cur==begin());
            container_.pop_front();
            return end();
          }
        else
          {
            container_.erase_after(prev);
            return prev;
          }
      }
  }


  inline void trig_slist::insert(trig_size_t index, mult_t multiplier)
  {
    find_result fr;
    find(index,fr);
    if (fr.second==end())
      {
        insert_new(index,multiplier,fr);
      }
    else
      {
        modify(multiplier,fr.first,fr.second);
      }
  }


  inline void trig_slist::invert_sign()
  {
    const iterator it_f=end();
    for (iterator it=begin();it!=it_f;++it)
      {
        it->second*=-1;
      }
  }


  // Probing implementations.
  inline double trig_slist::freq(const vector_psym_p &v) const
    {
      double retval=0.;
      const const_iterator it_f=end();
      for (const_iterator it=begin();it!=it_f;++it)
        {
          // We must be sure that there actually is a freq in every symbol we are going to use.
          if (v[it->first]->poly_eval().size()>1)
            {
              retval+=it->second*v[it->first]->poly_eval()[1];
            }
        }
      return retval;
    }


  inline double trig_slist::phase(const vector_psym_p &v) const
    {
      double retval=0.;
      const const_iterator it_f=end();
      for (const_iterator it=begin();it!=it_f;++it)
        {
          // We must be sure that there actually is a phase in every symbol we are going to use.
          if (v[it->first]->poly_eval().size()>0)
            {
              retval+=it->second*v[it->first]->poly_eval()[0];
            }
        }
      return retval;
    }


  inline double trig_slist::t_eval(const double &t, const vector_psym_p &v) const
    {
      double retval=0.;
      const const_iterator it_f=end();
      for (const_iterator it=begin();it!=it_f;++it)
        {
          retval+=it->second*v[it->first]->eval(t);
        }
      return retval;
    }


  inline short int trig_slist::sign() const
    {
      const const_iterator it_f=end();
      for (const_iterator it=begin();it!=it_f;++it)
        {
          if (it->second>0)
            {
              return 1;
            }
          if (it->second<0)
            {
              return -1;
            }
        }
      return 1;
    }


  inline bool trig_slist::operator>(const trig_slist &l2) const
    {
      const const_iterator it_f1=end(), it_f2=l2.end();
      const_iterator it1=begin(), it2=l2.begin();
      while (it1!=it_f1 && it2!=it_f2)
        {
          if (it1->first < it2->first)
            {
              return (it1->second > 0);
            }
          if (it2->first < it1->first)
            {
              return (it2->second < 0);
            }
          if (it1->second > it2->second)
            {
              return true;
            }
          if (it2->second > it1->second)
            {
              return false;
            }
          ++it1;
          ++it2;
        }
      if (it1!=it_f1)
        {
          return (it1->second > 0);
        }
      if (it2!=it_f2)
        {
          return (it2->second < 0);
        }
      return false;
    }


  inline bool trig_slist::operator<(const trig_slist &l2) const
    {
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


  inline size_t trig_slist::hasher() const
    {
      size_t seed=0;
      const const_iterator it_f=end();
      for (const_iterator it=begin();it!=it_f;++it)
        {
          boost::hash_combine(seed,it->first);
          boost::hash_combine(seed,it->second);
        }
      return seed;
    }


  inline bool trig_slist::is_zero() const
    {
      return empty();
    }


  inline bool trig_slist::smaller(const size_t &) const
    {
      return false;
    }


  inline bool trig_slist::larger(const size_t &) const
    {
      return false;
    }


  inline bool trig_slist::compatible(const size_t &) const
    {
      return true;
    }


  inline size_t trig_slist::data_footprint() const
    {
      return (sizeof(std::pair<trig_size_t,mult_t>)*container_.size());
    }


  inline bool trig_slist::checkup(const size_t &w) const
    {
      const const_iterator it_f=end();
      // Let's check there are not zero elements.
      for (const_iterator it=begin();it!=it_f;++it)
        {
          if (it->second==0)
            {
              std::cout << "Zero element in trig_slist." << std::endl;
              return false;
            }
        }
      // Now check that we do not have elements in excess.
      if (actual_width()>w)
        {
          std::cout << "Too many elements in trig_slist for given size." << std::endl;
          return false;
        }
      return true;
    }


  inline bool trig_slist::operator==(const trig_slist &l2) const
    {
      return (container_==l2.container_);
    }


  inline trig_slist &trig_slist::operator*=(mult_t n)
  {
    if (n==0)
      {
        container_.clear();
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


  template <short int Sign>
  inline trig_slist &trig_slist::algebraic_sum(const trig_slist &t2)
  {
    p_assert(this!=&t2);
    iterator it1=begin(), old_it1=end();
    const const_iterator it_f2=t2.end();
    const_iterator it2=t2.begin();
    while (it2!=it_f2)
      {
        // We are at the end of this, insert new term and do not increase iterators on this.
        if (it1==end())
          {
            old_it1=insert_after(it2->first,Sign*it2->second,old_it1);
            ++it2;
          }
        // We found a duplicate, modify it.
        else if (it1->first==it2->first)
          {
            old_it1=modify(it1->second+Sign*it2->second,old_it1,it1);
            // Take care in case we erased the element when modifying, we may have run out of elements.
            if (old_it1!=end())
              {
                it1=old_it1;
                ++it1;
              }
            else
              {
                it1=begin();
              }
            ++it2;
          }
        // We are past.
        else if (it1->first>it2->first)
          {
            old_it1=insert_after(it2->first,Sign*it2->second,old_it1);
            ++it2;
          }
        else if (it1->first<it2->first)
          {
            old_it1=it1;
            ++it1;
          }
      }
    return *this;
  }


  inline trig_slist &trig_slist::operator+=(const trig_slist &t2)
  {
    return algebraic_sum<+1>(t2);
  }


  inline trig_slist &trig_slist::operator-=(const trig_slist &t2)
  {
    return algebraic_sum<-1>(t2);
  }


  inline void trig_slist::trigmult(const trig_slist &l2, trig_slist &ret1, trig_slist &ret2) const
    {
      ret1=*this;
      ret1-=l2;
      ret2=*this;
      ret2+=l2;
    }


  // FIXME: maybe we can place this struct inside trig_slist and require that all trig_args provide
  // it, as an interface request? But then all classes that are expected to use the hasher should then
  // multi-inherit from a base hasher class? Otherwise how can we overload hash-value? Mh, maybe we
  // can specificy the hasher in multiindex containers. Investigate.
  /// Overload of hash_value function for piranha::trig_slist.
  /**
   * To be used in piranha::base_pseries for the hashed index.
   */
  inline size_t hash_value(const trig_slist &tm)
  {
    return tm.hasher();
  }
}

#endif // PIRANHA_TRIG_SLIST_H
