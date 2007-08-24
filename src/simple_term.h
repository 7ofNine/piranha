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

#ifndef PIRANHA_SIMPLE_TERM_H
#define PIRANHA_SIMPLE_TERM_H

#include "base_term.h"

namespace piranha
{
/// Simple Poisson series term class.
/**
 * Coefficient and trigonometric part are stored as statically allocated data members.
 */
  template <class Cf, class Trig>
    class simple_term:public base_term<Cf,Trig,simple_term<Cf,Trig> >
  {
    public:
/// Alias for self.
      typedef simple_term self;
/// Alias for ancestor.
      typedef base_term<Cf,Trig,simple_term<Cf,Trig> > ancestor;
/// Alias for coefficient type.
      typedef Cf cf_type;
/// Alias for trigonometric type.
      typedef Trig trig_type;
/// Default constructor.
      explicit simple_term():
      ancestor(true),private_cf_(),private_trig_()
        {}
// FIXME: replace bool with enum.
/// Constructor from coefficient and flavour.
      explicit simple_term(const cf_type &c, bool flavour=true):
      ancestor(flavour),private_cf_(c),private_trig_()
        {}
/// Generic builder.
/**
 * Build constructing coefficient from variable x, of type T.
 */
      template <class T>
        explicit simple_term(const T &x, bool flavour=true):
      ancestor(flavour),private_cf_(cf_type(x)),private_trig_()
        {}
/// Copy constructor from term with different coefficient type.
      template <class Cf2>
        explicit simple_term(const simple_term<Cf2,trig_type> &term):
      ancestor(term.g_flavour()),private_cf_(*term.g_cf()),private_trig_(*term.g_trig())
      {                                           /*BOOST_STATIC_ASSERT(sizeof(U)==0);*/
      }
// Getters
/// Get coefficient reference.
      cf_type *s_cf()
      {
        return &private_cf_;
      }
      const cf_type *g_cf() const
      {
        return &private_cf_;
      }
// TODO: rework this when flavour is included in trigs. Pay attention to const_mem_fun extractor
// in norm-based series indices.
/// Get flavour.
      bool &s_flavour()
      {
        return ancestor::private_flavour_;
      }
      const bool &g_flavour() const
      {
        return ancestor::private_flavour_;
      }
/// Get reference to trigonometric part.
      trig_type *s_trig()
      {
        return &private_trig_;
      }
      const trig_type *g_trig() const
      {
        return &private_trig_;
      }
      const trig_type &g_trig_ref() const
      {
        return private_trig_;
      }
      size_t footprint() const
      {
        return (sizeof(self)+g_trig()->data_footprint());
      }
/// Assignment operator.
      simple_term &operator=(const simple_term &t2)
      {
        if (this==&t2)
        {
          return *this;
        }
        *s_cf()=*t2.g_cf();
        *s_trig()=*t2.g_trig();
        ancestor::s_flavour()=t2.g_flavour();
        return *this;
      }
    private:
      cf_type     private_cf_;
      trig_type   private_trig_;
  };
}
#endif
