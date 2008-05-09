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

from _Core import *

import copy as __copy
import math as __math
try:
  import numpy as __numpy
  import matplotlib as __matplotlib
except:
  print "Warning: many capabilities of Pyranha rely on numpy and matplotlib. Please consider installing these packages."

def cos(arg):
  """
  Wrapper around standard cosine function. If arg provides a cos() method, it will be called, otherwise
  math.cos() will be called.
  """
  try:
    return arg.cos()
  except TypeError:
    return __math.cos(arg)

def sin(arg):
  """
  Wrapper around standard sine function. If arg provides a sin() method, it will be called, otherwise
  math.sin() will be called.
  """
  try:
    return arg.sin()
  except TypeError:
    return __math.sin(arg)

def besselJ(order,arg):
  """
  Wrapper around Core.besselJ. If arg provides a besselJ() method, it will be called, otherwise
  Core.besselJ will be called.
  """
  try:
    return arg.besselJ(order)
  except TypeError:
    return _Core.besselJ(order,arg)

def dbesselJ(order,arg):
  """
  Partial derivative of Bessel function of the first kind of integer order. It will call the dbesselJ() method
  of arg, if available.
  """
  try:
    return arg.dbesselJ(order)
  except AttributeError:
    raise AttributeError, "The dbesselJ() method is not available for this argument type, returning None."

def copy(arg):
  """Standard copy function. Lifted from the copy module."""
  return __copy.copy(arg)

def partial(arg,name):
  """
  Calculate partial derivative of arg with respect to argument name.

  Internally the partial() method of arg is called. If such method is not available, an AttributeError
  exception will be raised.
  """
  try:
    return arg.partial(name)
  except AttributeError:
    raise AttributeError, "The partial() method is not available for this argument type, returning None."

psym_manager = _Core.__psym_manager()
expo_truncator = _Core.__expo_truncator()
norm_truncator = _Core.__norm_truncator()

def tc(args, flambda, t0 = 0., t1 = None, step = None):
  try:
    args_list = list(args)
  except TypeError:
    args_list = [args]
  if (t1 == None and step != None) or (t1 != None and step == None):
    raise TypeError, 'You must specify both a final time AND a step size.'
  if t1 == None:
    args_list_eval = map(lambda x: x.eval(t0),args_list)
    return abs(flambda(*args_list).eval(t0)-flambda(*args_list_eval))
  else:
    if (t1 <= t0):
      raise ValueError, 't1 must be strictly greater than t0.'
    time_array = __numpy.arange(t0,t1,step)
    retval = __numpy.array([],dtype=float)
    retval.resize(len(time_array))
    result = flambda(*args_list);
    j=0
    for t in time_array:
      args_list_eval = map(lambda x: x.eval(t),args_list)
      retval[j]=abs(result.eval(t)-flambda(*args_list_eval))
      j=j+1
    return time_array,retval;
