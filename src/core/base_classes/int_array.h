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

#ifndef PIRANHA_INT_ARRAY_H
#define PIRANHA_INT_ARRAY_H

#include <boost/integer.hpp>
#include <boost/integer_traits.hpp>
#include <boost/functional/hash.hpp>
#include <boost/static_assert.hpp>
#include <boost/tuple/tuple.hpp>
#include <memory> // For std::allocator.

#include "../math.h" // For lg.
#include "../p_assert.h"

// Cast argument to int_array::max_fast_type pointer.
#define max_cast(arg) ((max_fast_type *)(arg))
#define derived_const_cast static_cast<Derived const *>(this)
#define derived_cast static_cast<Derived *>(this)
#define __PIRANHA_INT_ARRAY_TP_DECL int Bits, int Pos, bool Signed, class Allocator, class Derived
#define __PIRANHA_INT_ARRAY_TP Bits,Pos,Signed,Allocator,Derived

namespace piranha
{
  // Used below to define the integer type.
  template <int Bits, bool Signed>
    struct integer_type_chooser
  {
    BOOST_STATIC_ASSERT(Bits == 8 or Bits == 16);
    typedef typename boost::int_t<Bits>::fast type;
    typedef max_fast_int max_fast_type;
  };

  // Specialization for unsigned integers.
  template <int Bits>
    struct integer_type_chooser<Bits,false>
  {
    BOOST_STATIC_ASSERT(Bits == 8 or Bits == 16);
    typedef typename boost::uint_t<Bits>::fast type;
    typedef max_fast_uint max_fast_type;
  };

  /// Dynamically-sized integer array.
  /**
   * Parametrized to an integer sized Bits, which can be Signed or not. It contains also
   * a flavour boolean flag, which can be used in trigonometric parts
   * of Poisson series.
   */
  template <__PIRANHA_INT_ARRAY_TP_DECL>
    class int_array
  {
      typedef typename integer_type_chooser<Bits,Signed>::type value_type_;
      typedef typename integer_type_chooser<Bits,Signed>::max_fast_type max_fast_type;
      typedef typename Allocator::template rebind<value_type_>::other allocator_type;
    public:
      typedef value_type_ value_type;
      typedef uint8 size_type;
      BOOST_STATIC_ASSERT(sizeof(max_fast_type) % sizeof(value_type) == 0);
      BOOST_STATIC_ASSERT(Pos >= 0);
      static const int position = Pos;
      /// Default ctor.
      /**
       * Constructs empty array.
       */
      int_array():m_flavour(true),m_size(0),m_pack_size(0),m_ptr(allocator.allocate(0)) {}
      /// Copy ctor.
      int_array(const int_array &v):m_flavour(v.m_flavour),m_size(v.m_size),m_pack_size(v.m_pack_size),
        m_ptr(allocator.allocate(m_size))
      {
        packed_copy(m_ptr,v.m_ptr,m_size,m_pack_size);
      }
      /// Ctor from size.
      /**
       * Initialises to 0 all the elements.
       */
      int_array(const size_type &s):m_flavour(true),m_size(s),m_pack_size(s >> pack_shift),
        m_ptr(allocator.allocate(m_size))
      {
        size_type i;
        for (i=0;i < m_pack_size;++i)
        {
          max_cast(m_ptr)[i]=0;
        }
        for(i = i << pack_shift;i < m_size;++i)
        {
          m_ptr[i]=0;
        }
      }
      /// Dtor.
      ~int_array() {allocator.deallocate(m_ptr,m_size);}
      /// Print to stream
      void print(std::ostream &out_stream) const
      {
        stream_manager::setup_print(out_stream);
        for (size_t i=0;i < m_size;++i)
        {
          // We cast to max_fast_int, which is the largest integer type admitted..
          out_stream << (max_fast_int)(m_ptr[i]);
          // Print the separator iff this is not the last element.
          if (i != (size_t)(m_size-1))
          {
            out_stream << separator;
          }
        }
      }
      /// Array-like operator[], const version.
      const value_type &operator[](const size_t &n) const {return m_ptr[n];}
      /// Array-like operator[], mutable version.
      value_type &operator[](const size_t &n) {return m_ptr[n];}
      /// Return const reference to flavour.
      const bool &flavour() const {return m_flavour;}
      /// Return mutable reference to flavour.
      bool &flavour() {return m_flavour;}
      /// Return container size.
      size_t size() const {return m_size;}
      // TODO: should we pass size_t here, and test against size_type and throw in case of out-of-range request?
      /// Resize the container.
      /**
       * The existing elements are copied over.
       */
      void resize(const size_type &new_size)
      {
        value_type *new_ptr = alloc_if_size_differs(new_size);
        if (new_ptr == m_ptr)
        {
          return;
        }
        const size_type new_pack_size = (new_size >> pack_shift);
        // Copy to the minimum of the new sizes.
        packed_copy(new_ptr,m_ptr,std::min(m_size,new_size),std::min(m_pack_size,new_pack_size));
        // Zero the remaining elements, if any.
        for (size_type i=m_size;i < new_size;++i)
        {
          new_ptr[i]=0;
        }
        // Destroy old pointer and assign new data members.
        allocator.deallocate(m_ptr,m_size);
        m_ptr = new_ptr;
        m_size = new_size;
        m_pack_size = new_pack_size;
      }
      /// Assignment operator.
      int_array &operator=(const int_array &v)
      {
        value_type *new_ptr = alloc_if_size_differs(v.m_size);
        if (new_ptr != m_ptr)
        {
          allocator.deallocate(m_ptr,m_size);
          m_flavour = v.m_flavour;
          m_ptr = new_ptr;
          m_size = v.m_size;
          m_pack_size = v.m_pack_size;
        }
        packed_copy(m_ptr,v.m_ptr,m_size,m_pack_size);
        return *this;
      }
      /// Hash value.
      /**
       * Hashes only the integer elements of the array, not the flavour.
       */
      size_t hasher() const
      {
        size_t retval=0;
        size_type i;
        for (i=0;i < m_pack_size;++i)
        {
          boost::hash_combine(retval,max_cast(m_ptr)[i]);
        }
        for(i = i << pack_shift;i < m_size;++i)
        {
          boost::hash_combine(retval,m_ptr[i]);
        }
        return retval;
      }
      /// Equality test.
      /**
       * Tests only the integer elements of the array, not the flavour.
       */
      bool equal_to(const int_array &v) const
      {
        switch (m_size == v.m_size)
        {
          case true:
            size_type i;
            for (i=0;i < m_pack_size;++i)
            {
              if (max_cast(m_ptr)[i] != max_cast(v.m_ptr)[i])
              {
                return false;
              }
            }
            for(i = i << pack_shift;i < m_size;++i)
            {
              if (m_ptr[i] != v[i])
              {
                return false;
              }
            }
            return true;
          default:
            return false;
        }
      }
      bool lexicographic_less_than(const Derived &a2) const
      {
        p_assert(m_size == a2.m_size);
        for (size_t i=0; i < m_size; ++i)
        {
          if (m_ptr[i] < a2[i])
          {
            return true;
          }
          else if (m_ptr[i] > a2[i])
          {
            return false;
          }
        }
        return false;
      }
      /// Test for zero elements.
      /**
       * Returns true if all integer elements are zero, false otherwise.
       */
      bool is_zero() const
      {
        size_t i;
        for (i=0;i < m_pack_size;++i)
        {
          if (max_cast(m_ptr)[i] != 0)
          {
            return false;
          }
        }
        for(i = i << pack_shift;i < m_size;++i)
        {
          if (m_ptr[i] != 0)
          {
            return false;
          }
        }
        return true;
      }
      /// Do I need padding in order to be inserted in series?
      template <class ArgsTuple>
        bool needs_padding(const ArgsTuple &args_tuple) const
      {
        return (m_size < args_tuple.template get<Pos>().size());
      }
      /// Am I insertable in a series?
      template <class ArgsTuple>
        bool is_insertable(const ArgsTuple &args_tuple) const
      {
        return (m_size <= args_tuple.template get<Pos>().size());
      }
      template <class ArgsTuple>
        bool checkup(const ArgsTuple &args_tuple) const
      {
        switch (args_tuple.template get<Pos>().size() != m_size)
        {
          case true:
            std::cout << "Size mismatch in int_array." << std::endl;
            return false;
          default:
            return true;
        }
      }
      /// Pad right.
      template <class ArgsTuple>
        void pad_right(const ArgsTuple &args_tuple)
      {
        p_assert(args_tuple.template get<Pos>().size() >= m_size);
        resize(args_tuple.template get<Pos>().size());
      }
      template <class ArgsTuple, class Layout>
        void apply_layout(const ArgsTuple &, const Layout &l)
      {
        const size_t l_size = l.template get<Pos>().size();
        // The layout must have at least all arguments in this.
        p_assert(l_size >= m_size);
        // Memorize the old vector.
        const Derived old(*derived_const_cast);
        // Make space.
        resize(l_size);
        for (size_t i=0;i < l_size;++i)
        {
          switch (l.template get<Pos>()[i].first)
          {
            case true:
              p_assert(l.template get<Pos>()[i].second < old.m_size);
              m_ptr[i]=old[l.template get<Pos>()[i].second];
              break;
            case false:
              m_ptr[i]=0;
          }
        }
      }
    private:
      void packed_copy(value_type *new_ptr, const value_type *old_ptr, const size_type &size,
        const size_type &pack_size)
      {
        size_type i;
        for (i=0;i < pack_size;++i)
        {
          max_cast(new_ptr)[i]=max_cast(old_ptr)[i];
        }
        for(i = i << pack_shift;i < size;++i)
        {
          new_ptr[i]=old_ptr[i];
        }
      }
      // Allocate if new_size is different from current size.
      value_type *alloc_if_size_differs(const size_type &new_size)
      {
        switch (m_size == new_size)
        {
          case true:
            return m_ptr;
          default:
            return allocator.allocate(new_size);
        }
      }
    private:
      // Data members.
      /// Flavour.
      bool                    m_flavour;
      /// Size of the array.
      /**
       * Equal to the number of elements contained by the array.
       */
      size_type               m_size;
      /// Packed size of the array.
      /**
       * Defined by the integer division int_array::m_size / int_array::pack_mult.
       */
      size_type               m_pack_size;
      /// Pointer to the first value of the array.
      value_type              *m_ptr;
      /// Array allocator.
      static allocator_type   allocator;
      /// Pack multiplier.
      /**
       * Defined by the integer division between the number of bits of int_array::max_fast_type and the number
       * of bits of int_array::value_type.
       */
      static const size_type  pack_mult = sizeof(max_fast_type)/sizeof(value_type);
      /// Pack shifting.
      /**
       * Defined as the base-2 logarithm of int_array::pack_mult. If int_array::pack_mult is not a power of two,
       * a compilation error is produced.
       */
      static const size_type  pack_shift = lg<pack_mult>::value;
    public:
      static const char separator = ';';
  };

  template <__PIRANHA_INT_ARRAY_TP_DECL>
    typename int_array<__PIRANHA_INT_ARRAY_TP>::allocator_type int_array<__PIRANHA_INT_ARRAY_TP>::allocator;
};

#undef max_cast
#undef derived_const_cast
#undef derived_cast
#undef __PIRANHA_INT_ARRAY_TP_DECL
#undef __PIRANHA_INT_ARRAY_TP

#endif
