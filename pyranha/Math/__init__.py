# -*- coding: utf-8 -*-
# Copyright (C) 2007, 2008 by Francesco Biscani
# bluescarni@gmail.com
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the
# Free Software Foundation, Inc.,
# 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

from ._Math import *
from math import *   # that is a standard library

# NOTE: import cmath into cos, sin and allow for complex arguments?

def factorial(n):
    """
    Factorial of n.

    n must be a nonnegative (multiprecision) integer.
    """
    from ._Math import __factorial
    from pyranha.Core import integer
    if not isinstance(n,int) and not isinstance(n,integer): raise TypeError('n must be an integer')
    return __factorial(n)

def double_factorial(n):
    """
    Double factorial of n.

    n must be a nonnegative (multiprecision) integer.
    """
    from ._Math import __double_factorial
    from pyranha.Core import integer
    if not isinstance(n,int) and not isinstance(n,integer): raise TypeError('n must be an integer')
    return __double_factorial(n)

def r_factorial(x,n):
    """
    Rising factorial x^(n).

    x must be a numerical type, n must be a nonnegative integer.
    """
    from ._Math import __r_factorial
    if not isinstance(n,int): raise TypeError('n must be an integer')
    return __r_factorial(x,n)

def f_factorial(x,n):
    """
    Falling factorial (x)_n.

    x must be a numerical type, n must be a nonnegative integer.
    """
    from ._Math import __f_factorial
    if not isinstance(n,int): raise TypeError('n must be an integer')
    return __f_factorial(x,n)

def choose(x,k):
    """
    Generalised binomial coefficient (x over k).

    k must be an integer, while x can be any (complex) numerical type. This implementation returns
    0 whenever k < 0 or k > n.
    """
    from ._Math import __choose
    if not isinstance(k,int): raise TypeError('k must be an int')
    return __choose(x,k)

def ei(arg):
    """
    Complex exponential of arg.
    """
    import math
    try:
        return arg.ei()
    except AttributeError:
        return complex(math.cos(arg),math.sin(arg))

def cos(arg):
    """
    Wrapper around standard cosine function. If arg provides a cos() method, it will be called, otherwise
    math.cos() will be called.
    """
    import math
    try:
        return arg.cos()
    except AttributeError:
        return math.cos(arg)

def sin(arg):
    """
    Wrapper around standard sine function. If arg provides a sin() method, it will be called, otherwise
    math.sin() will be called.
    """
    import math
    try:
        return arg.sin()
    except AttributeError:
        return math.sin(arg)

def root(n,arg):
    """
    n-th root of argument arg. If arg provides a root() method, it will be called, otherwise
    the standard ** operator will be called.
    """
    if not isinstance(n,int): raise TypeError('n must be an int')
    try:
        return arg.root(n)
    except AttributeError:
        return arg**(1./n)

def hyperF(a_sequence,b_sequence,z,limit = None):
    """
    Hypergeometric series of argument z, with the a and b parameters provided as sequences
    of objects that can be used to construct a rational.
    """
    from pyranha.Core import rational
    if not limit is None and not isinstance(limit,int):
        raise TypeError('limit must be either None or an integer value.')
    try:
        if limit is None:
            return z.hyperF([rational(a) for a in a_sequence],[rational(b) for b in b_sequence])
        else:
            return z.hyperF([rational(a) for a in a_sequence],[rational(b) for b in b_sequence],limit)
    except TypeError as ArgumentError:
        raise TypeError('inputs a_sequence and b_sequence must be sequences of elements from which rationals can be constructed.')
    except AttributeError:
        raise TypeError('z does not provide an hyperF() method.')

def dhyperF(a_sequence,b_sequence,z,limit = None,order = 1):
    """
    Derivative of the hypergeometric series of argument z, with the a and b parameters provided as sequences
    of objects that can be used to construct a rational. Order of derivation defaults to one.
    """
    from pyranha.Core import rational
    if not limit is None and not isinstance(limit,int):
        raise TypeError('limit must be either None or an integer value.')
    if not isinstance(order,int):
        raise TypeError('order must be an integer value.')
    try:
        if limit is None:
            return z.dhyperF(order,[rational(a) for a in a_sequence],[rational(b) for b in b_sequence])
        else:
            return z.dhyperF(order,[rational(a) for a in a_sequence],[rational(b) for b in b_sequence],limit)
    except TypeError as ArgumentError:
        raise TypeError('inputs a_sequence and b_sequence must be sequences of elements from which rationals can be constructed.')
    except AttributeError:
        raise TypeError('z does not provide an hyperF() method.')

def besselJ(order,arg):
    """
    Bessel function of the first kind of integer order of argument arg.
    If arg provides a besselJ() method, it will be called, otherwise _Math.besselJ will be called.
    """
    if not isinstance(order,int): raise TypeError('order must be an int')
    try:
        return arg.besselJ(order)
    except AttributeError:
        return _Math.besselJ(order,arg)

def dbesselJ(order,arg):
    """
    Partial derivative of Bessel function of the first kind of integer order of argument arg.
    It will call the dbesselJ() method of arg, if available, otherwise an AttributeError exception
    will be raised.
    """
    if not isinstance(order,int): raise TypeError('order must be an int')
    try:
        return arg.dbesselJ(order)
    except AttributeError:
        return _Math.dbesselJ(order,arg)

def besselJ_div_m(order,arg,m = 1):
    """
    Bessel function of the first kind of integer order of argument arg divided by arg**m.
    If arg provides a besselJ_div_m() method, it will be called, otherwise _Math.besselJ will be called.
    """
    if not isinstance(order,int) or not isinstance(m,int): raise TypeError('order and m must be ints')
    try:
        return arg.besselJ_div_m(order,m)
    except AttributeError:
        return _Math.besselJ(order,arg) / arg**m

def legendrePnm(n,m,arg,*extra_arg):
    """
    Associated Legendre function of integer degree n and order m of argument arg. The optional argument extra_arg
    can be used to provide a value for sqrt(1-arg**2), which will be used for the calculation of
    Pnm using recurrence relations. If extra_arg is not provided, the function will try to calculate
    sqrt(1-arg**2) on its own.
    """
    if not isinstance(n,int) or not isinstance(m,int): raise TypeError('n and m must be ints')
    if len(list(extra_arg)) > 1:
        raise TypeError("Please provide at most one extra argument.")
    try:
        return arg.legendrePnm(n,m,*extra_arg)
    except AttributeError:
        return _Math.legendrePnm(n,m,arg)

def legendrePn(n,arg):
    """
    Legendre polynomial of integer degree n of argument arg.
    """
    if not isinstance(n,int): raise TypeError('n must be an int')
    try:
        return arg.legendrePn(n)
    except AttributeError:
        return _Math.legendrePn(n,arg)

def Ynm(n, m, theta, phi, emi_phi = None, alpha = None, beta = None, gamma = None):
    """
    Non-normalised spherical harmonic of integer degree n and order m of the colatitude theta and longitude phi.
    
    If at least one of alpha, beta or gamma is provided, the spherical harmonic will be rotated under
    the 3-1-3 Euler angles alpha, beta and gamma using Wigner's theorem for the rotation of spherical harmonics.
    
    If emi_phi is provided, it will be used internally as exp(-i*phi) and phi will be assumed to mean exp(i*phi),
    otherwise the routine will assume that phi is the longitude and it will try to calculate exp(-i*phi) on its own.

    Other than n and m, which must always be integers, the arguments of this function must be of homogeneous type
    (i.e., all series or all numbers).
    """
    if not isinstance(n,int) or not isinstance(m,int): raise TypeError('n and m must be ints')
    if alpha != None or beta != None or gamma != None:
        a = alpha or type(theta)(0)
        b = beta or type(theta)(0)
        g = gamma or type(theta)(0)
        if emi_phi == None:
            try:
                return type(theta).Ynm(n,m,theta,phi,a,b,g)
            except AttributeError:
                return _Math.Ynm(n,m,theta,phi,a,b,g)
        else:
            ei_phi = phi
            try:
                return type(theta).Ynm(n,m,theta,ei_phi,emi_phi,a,b,g)
            except AttributeError:
                return _Math.Ynm(n,m,theta,ei_phi,emi_phi,a,b,g)
    else:
        if emi_phi == None:
            try:
                return type(theta).Ynm(n,m,theta,phi)
            except AttributeError:
                return _Math.Ynm(n,m,theta,phi)
        else:
            ei_phi = phi
            try:
                return type(theta).Ynm(n,m,theta,ei_phi,emi_phi)
            except AttributeError:
                return _Math.Ynm(n,m,theta,ei_phi,emi_phi)

def partial(arg,p,n = 1):
    """
    Calculate the n-th partial derivative of arg with respect to symbol named p.

    Internally the partial() method of arg is called. If such method is not available, an AttributeError
    exception will be raised.
    """
    if not isinstance(n,int): raise TypeError('n must be an int')
    try:
        return arg.partial(p,n)
    except AttributeError:
        raise AttributeError("The partial() method is not available for argument type '" + str(type(arg)) + "'")

def integrate(arg,p):
    """
    Indefinite integral of arg with respect to symbol named p.

    Internally the integrate() method of arg is called. If such method is not available, an AttributeError
    exception will be raised.
    """
    try:
        return arg.integrate(p)
    except AttributeError:
        raise AttributeError("The integrate() method is not available for argument type '" + str(type(arg)) + "'")

def einpi2(n):
    """
    Complex exponential of n*pi/2.
    """
    if not isinstance(n,int): raise TypeError('n must be an int')
    if n & 1:
        if (n - 1) & 3:
            return complex(0,-1)
        else:
            return complex(0,1)
    else:
        if n & 3:
            return complex(-1,0)
        else:
            return complex(1,0)

def grad(l_var,s):
    """
    Gradient of series s with respect to variables whose names are listed in l_var.
    """
    for n in l_var:
        if not isinstance(n,str):
            raise TypeError('Please pass a list of variable names as first argument.')
    return [partial(s,n) for n in l_var]
