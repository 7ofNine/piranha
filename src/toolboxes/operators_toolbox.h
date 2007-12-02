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

#ifndef PIRANHA_OPERATORS_TOOLBOX_H
#define PIRANHA_OPERATORS_TOOLBOX_H

namespace piranha
{
  template <class Derived>
    class real_operators_toolbox
  {
    public:
      Derived &operator=(const Derived &p)
      {
        return static_cast<Derived *>(this)->assign(p);
      }
      Derived &operator+=(const double &x)
      {
        return static_cast<Derived *>(this)->add_generic(x);
      }
      Derived &operator+=(const Derived &p)
      {
        return static_cast<Derived *>(this)->add(p);
      }
      Derived operator+(const double &x) const
      {
        Derived retval(*static_cast<Derived const *>(this));
        retval+=x;
        return retval;
      }
      Derived operator+(const Derived &p) const
      {
        Derived retval(*static_cast<Derived const *>(this));
        retval+=p;
        return retval;
      }
      Derived &operator-=(const double &x)
      {
        return static_cast<Derived *>(this)->add_generic(-x);
      }
      Derived &operator-=(const Derived &p)
      {
        return static_cast<Derived *>(this)->subtract(p);
      }
      Derived operator-(const double &x) const
      {
        Derived retval(*static_cast<Derived const *>(this));
        retval-=x;
        return retval;
      }
      Derived operator-(const Derived &p) const
      {
        Derived retval(*static_cast<Derived const *>(this));
        retval-=p;
        return retval;
      }
      Derived &operator*=(int n)
      {
        return static_cast<Derived *>(this)->mult_by(n);
      }
      Derived &operator*=(const Derived &p)
      {
        return static_cast<Derived *>(this)->mult_by_self(p);
      }
      Derived &operator*=(const double &x)
      {
        return static_cast<Derived *>(this)->mult_by(x);
      }
      Derived operator*(int n) const
      {
        Derived retval(*static_cast<Derived const *>(this));
        retval*=n;
        return retval;
      }
      Derived operator*(const Derived &p) const
      {
        Derived retval(*static_cast<Derived const *>(this));
        retval*=p;
        return retval;
      }
      Derived &operator/=(int n)
      {
        return static_cast<Derived *>(this)->divide_by(n);
      }
      Derived &operator/=(const double &x)
      {
        return static_cast<Derived *>(this)->divide_by(x);
      }
      Derived operator/(const double &x) const
      {
        Derived retval(*static_cast<Derived const *>(this));
        retval/=x;
        return retval;
      }
      Derived operator/(int n) const
      {
        Derived retval(*static_cast<Derived const *>(this));
        retval/=n;
        return retval;
      }
  };

  template <class real_Derived>
    class complex_operators_toolbox
  {
    public:
      typedef std::complex<real_Derived> Derived;
      Derived &operator=(const Derived &p)
      {
        return static_cast<Derived *>(this)->assign(p);
      }
      Derived &operator+=(const double &x)
      {
        return static_cast<Derived *>(this)->add_generic(x);
      }
      Derived &operator+=(const Derived &p)
      {
        return static_cast<Derived *>(this)->add(p);
      }
      Derived operator+(const double &x) const
      {
        Derived retval(*static_cast<Derived const *>(this));
        retval+=x;
        return retval;
      }
      Derived operator+(const Derived &p) const
      {
        Derived retval(*static_cast<Derived const *>(this));
        retval+=p;
        return retval;
      }
      Derived &operator-=(const double &x)
      {
        return static_cast<Derived *>(this)->add_generic(-x);
      }
      Derived &operator-=(const Derived &p)
      {
        return static_cast<Derived *>(this)->subtract(p);
      }
      Derived operator-(const double &x) const
      {
        Derived retval(*static_cast<Derived const *>(this));
        retval-=x;
        return retval;
      }
      Derived operator-(const Derived &p) const
      {
        Derived retval(*static_cast<Derived const *>(this));
        retval-=p;
        return retval;
      }
      Derived &operator*=(int n)
      {
        static_cast<Derived *>(this)->mult_by_int(n);
        return *static_cast<Derived *>(this);
      }
      Derived &operator*=(const Derived &p)
      {
        return static_cast<Derived *>(this)->mult_by_self(p);
      }
      Derived &operator*=(const double &x)
      {
        static_cast<Derived *>(this)->mult_by_double(x);
        return *static_cast<Derived *>(this);
      }
      Derived operator*(int n) const
      {
        Derived retval(*static_cast<Derived const *>(this));
        retval*=n;
        return retval;
      }
      Derived operator*(const Derived &p) const
      {
        Derived retval(*static_cast<Derived const *>(this));
        retval*=p;
        return retval;
      }
      Derived &operator/=(int n)
      {
        static_cast<Derived *>(this)->generic_division(n);
        return *static_cast<Derived *>(this);
      }
      Derived &operator/=(const double &x)
      {
        return *static_cast<Derived *>(this)*=(1./x);
      }
      Derived operator/(const double &x) const
      {
        Derived retval(*static_cast<Derived const *>(this));
        retval/=x;
        return retval;
      }
      Derived operator/(int n) const
      {
        Derived retval(*static_cast<Derived const *>(this));
        retval/=n;
        return retval;
      }
// Complex specifics.
      Derived &operator=(const real_Derived &p)
      {
        return static_cast<Derived *>(this)->assign(p);
      }
      Derived &operator+=(const real_Derived &p)
      {
        return static_cast<Derived *>(this)->add(p);
      }
      Derived operator+(const real_Derived &p) const
      {
        Derived retval(*static_cast<Derived const *>(this));
        retval+=p;
        return retval;
      }
      Derived &operator-=(const real_Derived &p)
      {
        return static_cast<Derived *>(this)->add(p,false);
      }
      Derived operator-(const real_Derived &p) const
      {
        Derived retval(*static_cast<Derived const *>(this));
        retval-=p;
        return retval;
      }
      Derived &operator*=(const real_Derived &p)
      {
        return static_cast<Derived *>(this)->mult_by_self(p);
      }
      Derived operator*(const real_Derived &p) const
      {
        Derived retval(*static_cast<Derived const *>(this));
        retval*=p;
        return retval;
      }
  };
}
#endif
