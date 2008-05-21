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

from _Math import *

import pyranha as __pyranha
import math as __math

def complex(arg):
  """
  Wrapper around standard builtin "complex". If arg provides a complex() method, it will be called, otherwise
  builtin complex will be called.
  """
  try:
    return arg.complex()
  except AttributeError:
    return complex(arg)

def cos(arg):
  """
  Wrapper around standard cosine function. If arg provides a cos() method, it will be called, otherwise
  math.cos() will be called.
  """
  try:
    return arg.cos()
  except AttributeError:
    return __math.cos(arg)

def sin(arg):
  """
  Wrapper around standard sine function. If arg provides a sin() method, it will be called, otherwise
  math.sin() will be called.
  """
  try:
    return arg.sin()
  except AttributeError:
    return __math.sin(arg)

def root(n,arg):
  """
  Wrapper around root. If arg provides a root() method, it will be called, otherwise
  the standard ** operator will be called.
  """
  try:
    return arg.root(n)
  except AttributeError:
    return arg**(1./n)

def besselJ(order,arg):
  """
  Wrapper around _Math.besselJ. If arg provides a besselJ() method, it will be called, otherwise
  _Math.besselJ will be called.
  """
  try:
    return arg.besselJ(order)
  except AttributeError:
    return _Math.besselJ(order,arg)

def dbesselJ(order,arg):
  """
  Partial derivative of Bessel function of the first kind of integer order. It will call the dbesselJ() method
  of arg, if available.
  """
  try:
    return arg.dbesselJ(order)
  except AttributeError:
    raise AttributeError, "The dbesselJ() method is not available for this argument type, returning None."

def partial(arg,name):
  """
  Calculate partial derivative of arg with respect to argument name.

  Internally the partial() method of arg is called. If such method is not available, an AttributeError
  exception will be raised.
  """
  try:
    return arg.partial(name)
  except AttributeError:
    raise AttributeError, "The partial() method is not available for this argument type."