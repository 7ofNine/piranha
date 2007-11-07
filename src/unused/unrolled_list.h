#ifndef PIRANHA_UNROLLED_LIST
#define PIRANHA_UNROLLED_LIST

#include <boost/static_assert.hpp>
#include <cmath>
#include <cstring>
#include <ext/array_allocator.h>
#include <ext/bitmap_allocator.h>
#include <ext/mt_allocator.h>
#include <ext/pool_allocator.h>
#include <limits.h>

namespace piranha
{
/// Unrolled singly linked list.
  template <class T, int UnrollN>
    class unrolled_list
  {
    private:
      typedef unsigned short int num_t;
// Check that unroll limit is sane.
      BOOST_STATIC_ASSERT(UnrollN >= 2);
      BOOST_STATIC_ASSERT(UnrollN <= USHRT_MAX/2);
// Check that UnrollN is a multiple of 2.
      BOOST_STATIC_ASSERT((UnrollN/2)*2 == UnrollN);
      class node
      {
        private:
          friend class unrolled_list;
          friend class iterator;
          explicit node()
          {BOOST_STATIC_ASSERT(sizeof(T) == 0);}
          void assign_new(const T &element)
          {
            next=0;
            num_elements=1;
            elements[0]=element;
          }
          void push_front(const T &element)
          {
// TODO: assert num_elements<=UnrollN and num_elements>0.
            if (num_elements < UnrollN)
            {
              memmove((T *)elements+1,(T const *)elements,sizeof(T)*num_elements);
              elements[0]=element;
              ++num_elements;
            }
            else
            {
// Create new node and setup its properties.
              node *new_node=pool_node_.allocate(1);
              new_node->num_elements=UnrollN/2;
              new_node->next=next;
// Copy second half of current elements to the first half of new node.
              memcpy((T *)new_node->elements,(T const *)elements+UnrollN/2,sizeof(T)*UnrollN/2);
// Move up by one current elements that have to be retained (i.e., first half).
              memmove((T *)elements+1,(T const *)elements,sizeof(T)*UnrollN/2);
// Prepend the desired element.
              elements[0]=element;
// Current elements number has to be updated.
              num_elements=UnrollN/2+1;
// Connect to new node.
              next=new_node;
            }
          }
        private:
          num_t num_elements;
          node  *next;
          T     elements[UnrollN];
      };
      typedef __gnu_cxx::__mt_alloc<node>  /*__gnu_cxx::__pool_alloc<node>*/ pool_node;
    public:
      class iterator
      {
        public:
          explicit iterator(node *n, num_t m):private_position_(m),private_node_(n)
          {}
          iterator &operator++()
          {
            ++private_position_;
            if (private_position_ == private_node_->num_elements)
            {
              private_node_ = private_node_->next;
              private_position_ = 0;
            }
            return *this;
          }
          bool operator!=(const iterator &it) const
          {
            return ((private_node_ != it.private_node_) || (private_position_ != it.private_position_));
          }
          const T &operator*() const
          {
            return (private_node_->elements[private_position_]);
          }
          T &operator*()
          {
            return (private_node_->elements[private_position_]);
          }
        private:
          explicit iterator()
          {}
        private:
          num_t private_position_;
          node  *private_node_;
      };
    public:
      unrolled_list():private_head_(0)
      {}
      ~unrolled_list()
      {
        if (!empty())
        {
          node *tmp;
          while (private_head_->next != 0)
          {
            tmp=private_head_->next;
            pool_node_.deallocate(private_head_,1);
            private_head_=tmp;
          }
          pool_node_.deallocate(private_head_,1);
        }
      }
      iterator begin()
      {
        return iterator(private_head_,0);
      }
      iterator end()
      {
        return iterator(0,0);
      }
      void push_front(const T &element)
      {
        if (empty())
        {
          private_head_=pool_node_.allocate(1);
          private_head_->assign_new(element);
        }
        else
        {
          private_head_->push_front(element);
        }
      }
      bool empty() const
      {
        return (private_head_ == 0);
      }
    private:
      node                            *private_head_;
      static pool_node                pool_node_;
  };

  template <class T, int UnrollN>
    typename unrolled_list<T,UnrollN>::pool_node unrolled_list<T,UnrollN>::pool_node_;
}

#endif
