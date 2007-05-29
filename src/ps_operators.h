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

#ifndef PIRANHA_PS_OPERATORS_H
#define PIRANHA_PS_OPERATORS_H

// Macro for the definition of operators in pseries declarations
#define __PS_OPERATORS(ps_type,ancestor_type) \
    ps_type &operator=(const ps_type &p)\
    {\
        ancestor_type::basic_assignment(p);\
        return *this;\
    }\
    template <class Cf2,class Trig2,template <class,class> class I2>\
        ps_type &operator+=(const ps_type<Cf2,Trig2,I2> &p)\
    {\
        ancestor_type::merge_with(p);\
        return *this;\
    }\
    template <class Cf2,class Trig2,template <class,class> class I2>\
        ps_type &operator-=(const ps_type<Cf2,Trig2,I2> &p)\
    {\
        ancestor_type::merge_with(p,false);\
        return *this;\
    }\
    template <class Cf2,class Trig2,template <class,class> class I2>\
        ps_type &operator*=(const ps_type<Cf2,Trig2,I2> &p)\
    {\
        ancestor_type::basic_ps_mult(p);\
        return *this;\
    }\
    template <class Cf2,class Trig2,template <class,class> class I2>\
        ps_type &operator/=(const ps_type<Cf2,Trig2,I2> &p)\
    {\
        *this*=(p.pow(-1));\
        return *this;\
    }\
    template <class Cf2,class Trig2,template <class,class> class I2>\
        ps_type operator+(const ps_type<Cf2,Trig2,I2> &p) const\
    {\
        ps_type retval(*this);\
        retval+=p;\
        return retval;\
    }\
    template <class Cf2,class Trig2,template <class,class> class I2>\
        ps_type operator-(const ps_type<Cf2,Trig2,I2> &p) const\
    {\
        ps_type retval(*this);\
        retval-=p;\
        return retval;\
    }\
    template <class Cf2,class Trig2,template <class,class> class I2>\
        ps_type operator*(const ps_type<Cf2,Trig2,I2> &p) const\
    {\
        ps_type retval(*this);\
        retval*=p;\
        return retval;\
    }\
    template <class Cf2,class Trig2,template <class,class> class I2>\
        ps_type operator/(const ps_type<Cf2,Trig2,I2> &p) const\
    {\
        ps_type retval(*this);\
        retval/=p;\
        return retval;\
    }\
    ps_type &operator+=(int x)\
    {\
        ancestor_type::generic_merge(x);\
        return *this;\
    }\
    ps_type &operator+=(const double &x)\
    {\
        ancestor_type::generic_merge(x);\
        return *this;\
    }\
    ps_type &operator+=(const long double &x)\
    {\
        ancestor_type::generic_merge(x);\
        return *this;\
    }\
    ps_type operator+(const double &x) const\
    {\
        ps_type retval(*this);\
        retval+=x;\
        return retval;\
    }\
    ps_type &operator-=(int x)\
    {\
        ancestor_type::generic_merge(-x);\
        return *this;\
    }\
    ps_type &operator-=(const double &x)\
    {\
        ancestor_type::generic_merge(-x);\
        return *this;\
    }\
    ps_type &operator-=(const long double &x)\
    {\
        ancestor_type::generic_merge(-x);\
        return *this;\
    }\
    ps_type operator-(const double &x) const\
    {\
        ps_type retval(*this);\
        retval-=x;\
        return retval;\
    }\
    ps_type &operator*=(int n)\
    {\
        ancestor_type::mult_by_int(n);\
        return *this;\
    }\
    ps_type &operator*=(const double &x)\
    {\
        ancestor_type::generic_mult(x);\
        return *this;\
    }\
    ps_type &operator*=(const long double &x)\
    {\
        ancestor_type::generic_mult(x);\
        return *this;\
    }\
    ps_type operator*(const double &x) const\
    {\
        ps_type retval(*this);\
        retval*=x;\
        return retval;\
    }\
    ps_type operator*(int n) const\
    {\
        ps_type retval(*this);\
        retval*=n;\
        return retval;\
    }\
    ps_type &operator/=(const double &x)\
    {\
        return *this*=(1./x);\
    }\
    ps_type &operator/=(const long double &x)\
    {\
        return *this*=(1./x);\
    }\
    ps_type operator/(const double &x) const\
    {\
        ps_type retval(*this);\
        retval/=x;\
        return retval;\
    }\
    ps_type pow(const double &power) const\
    {\
        ps_type retval(*this);\
        retval.basic_pow(power);\
        return retval;\
    }

#endif // PIRANHA_PS_OPERATORS_H
