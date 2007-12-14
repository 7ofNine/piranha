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
#define __PIRANHA_BASE_OPERATORS Derived &operator=(const Derived &p)\
{\
  return static_cast<Derived *>(this)->assign(p);\
}\
Derived &operator=(int n)\
{\
  return static_cast<Derived *>(this)->assign(n);\
}\
Derived &operator=(const double &x)\
{\
  return static_cast<Derived *>(this)->assign(x);\
}\
Derived &operator+=(int n)\
{\
  return static_cast<Derived *>(this)->add(n);\
}\
Derived &operator+=(const double &x)\
{\
  return static_cast<Derived *>(this)->add(x);\
}\
Derived &operator+=(const Derived &p)\
{\
  return static_cast<Derived *>(this)->add(p);\
}\
Derived operator+(int n) const\
{\
  Derived retval(*static_cast<Derived const *>(this));\
  retval+=n;\
  return retval;\
}\
Derived operator+(const double &x) const\
{\
  Derived retval(*static_cast<Derived const *>(this));\
  retval+=x;\
  return retval;\
}\
Derived operator+(const Derived &p) const\
{\
  Derived retval(*static_cast<Derived const *>(this));\
  retval+=p;\
  return retval;\
}\
Derived &operator-=(int n)\
{\
  return static_cast<Derived *>(this)->subtract(n);\
}\
Derived &operator-=(const double &x)\
{\
  return static_cast<Derived *>(this)->subtract(x);\
}\
Derived &operator-=(const Derived &p)\
{\
  return static_cast<Derived *>(this)->subtract(p);\
}\
Derived operator-(int n) const\
{\
  Derived retval(*static_cast<Derived const *>(this));\
  retval-=n;\
  return retval;\
}\
Derived operator-(const double &x) const\
{\
  Derived retval(*static_cast<Derived const *>(this));\
  retval-=x;\
  return retval;\
}\
Derived operator-(const Derived &p) const\
{\
  Derived retval(*static_cast<Derived const *>(this));\
  retval-=p;\
  return retval;\
}\
Derived &operator*=(int n)\
{\
  return static_cast<Derived *>(this)->mult_by(n);\
}\
Derived &operator*=(const double &x)\
{\
  return static_cast<Derived *>(this)->mult_by(x);\
}\
Derived &operator*=(const Derived &p)\
{\
  return static_cast<Derived *>(this)->mult_by(p);\
}\
Derived operator*(int n) const\
{\
  Derived retval(*static_cast<Derived const *>(this));\
  retval*=n;\
  return retval;\
}\
Derived operator*(const double &x) const\
{\
  Derived retval(*static_cast<Derived const *>(this));\
  retval*=x;\
  return retval;\
}\
Derived operator*(const Derived &p) const\
{\
  Derived retval(*static_cast<Derived const *>(this));\
  retval*=p;\
  return retval;\
}\
Derived &operator/=(int n)\
{\
  return static_cast<Derived *>(this)->divide_by(n);\
}\
Derived &operator/=(const double &x)\
{\
  return static_cast<Derived *>(this)->divide_by(x);\
}\
Derived operator/(int n) const\
{\
  Derived retval(*static_cast<Derived const *>(this));\
  retval/=n;\
  return retval;\
}\
Derived operator/(const double &x) const\
{\
  Derived retval(*static_cast<Derived const *>(this));\
  retval/=x;\
  return retval;\
}

  template <class Derived>
    class real_operators_toolbox
  {
    public:
      __PIRANHA_BASE_OPERATORS
  };

// Requires: complex toolbox.
  template <class real_Derived>
    class complex_operators_toolbox
  {
      typedef std::complex<real_Derived> Derived;
    public:
      __PIRANHA_BASE_OPERATORS;
// Complex specifics.
      Derived &operator=(const std::complex<int> &n)
      {
        return static_cast<Derived *>(this)->complex_assign(n);
      }
      Derived &operator=(const std::complex<double> &x)
      {
        return static_cast<Derived *>(this)->complex_assign(x);
      }
      Derived &operator=(const real_Derived &p)
      {
        return static_cast<Derived *>(this)->complex_assign(p);
      }
      Derived &operator+=(const std::complex<int> &n)
      {
        return static_cast<Derived *>(this)->complex_add(n);
      }
      Derived &operator+=(const std::complex<double> &x)
      {
        return static_cast<Derived *>(this)->complex_add(x);
      }
      Derived &operator+=(const real_Derived &p)
      {
        return static_cast<Derived *>(this)->complex_add(p);
      }
      Derived operator+(const std::complex<int> &n) const
      {
        Derived retval(*static_cast<Derived const *>(this));
        retval+=n;
        return retval;
      }
      Derived operator+(const std::complex<double> &x) const
      {
        Derived retval(*static_cast<Derived const *>(this));
        retval+=x;
        return retval;
      }
      Derived operator+(const real_Derived &p) const
      {
        Derived retval(*static_cast<Derived const *>(this));
        retval+=p;
        return retval;
      }
      Derived &operator-=(const std::complex<int> &n)
      {
        return static_cast<Derived *>(this)->complex_subtract(n);
      }
      Derived &operator-=(const std::complex<double> &x)
      {
        return static_cast<Derived *>(this)->complex_subtract(x);
      }
      Derived &operator-=(const real_Derived &p)
      {
        return static_cast<Derived *>(this)->complex_subtract(p);
      }
      Derived operator-(const std::complex<int> &n) const
      {
        Derived retval(*static_cast<Derived const *>(this));
        retval-=n;
        return retval;
      }
      Derived operator-(const std::complex<double> &x) const
      {
        Derived retval(*static_cast<Derived const *>(this));
        retval-=x;
        return retval;
      }
      Derived operator-(const real_Derived &p) const
      {
        Derived retval(*static_cast<Derived const *>(this));
        retval-=p;
        return retval;
      }
      Derived &operator*=(const std::complex<int> &n)
      {
        return static_cast<Derived *>(this)->complex_mult_by(n);
      }
      Derived &operator*=(const std::complex<double> &x)
      {
        return static_cast<Derived *>(this)->complex_mult_by(x);
      }
      Derived &operator*=(const real_Derived &p)
      {
        return static_cast<Derived *>(this)->complex_mult_by(p);
      }
      Derived operator*(const std::complex<int> &n) const
      {
        Derived retval(*static_cast<Derived const *>(this));
        retval*=n;
        return retval;
      }
      Derived operator*(const std::complex<double> &x) const
      {
        Derived retval(*static_cast<Derived const *>(this));
        retval*=x;
        return retval;
      }
      Derived operator*(const real_Derived &p) const
      {
        Derived retval(*static_cast<Derived const *>(this));
        retval*=p;
        return retval;
      }
      Derived &operator/=(const std::complex<int> &n)
      {
        return static_cast<Derived *>(this)->complex_divide_by(n);
      }
      Derived &operator/=(const std::complex<double> &x)
      {
        return static_cast<Derived *>(this)->complex_divide_by(x);
      }
      Derived operator/(const std::complex<int> &n) const
      {
        Derived retval(*static_cast<Derived const *>(this));
        retval/=n;
        return retval;
      }
      Derived operator/(const std::complex<double> &x) const
      {
        Derived retval(*static_cast<Derived const *>(this));
        retval/=x;
        return retval;
      }
  };
#undef __PIRANHA_BASE_OPERATORS
}
#endif
