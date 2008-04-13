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

from @MODULE_NAME@ import *

import copy as _copy
import math as _math
import re as _re

try:
  import pylab as _pylab
  def __plot_range_evaluator(self,*args,**kwargs):
    _pylab.plot(self.times(),self.values(),*args,**kwargs)
  def __plot_tc(self,*args,**kwargs):
    values = list()
    for i in range(self.plain_eval().times().__len__()):
      values.append(abs(self.plain_eval().values()[i]-self.manip_eval().values()[i]))
    _pylab.semilogy(self.plain_eval().times(),self.plain_eval().values(),*args,**kwargs)
    _pylab.semilogy(self.manip_eval().times(),values,*args,**kwargs)
  range_evaluator.plot = __plot_range_evaluator
  def __add_plot_to_tcs():
    tc_list=[i for i in dir(@MODULE_NAME@) if _re.search('^tc_.*',i)]
    for i in tc_list:
      eval(i).plot = __plot_tc
  __add_plot_to_tcs()
except ImportError:
  print "Matplotlib is not installed, disabling plotting methods."
except NameError:
  pass

# Handy definitions of common mathematical functions: try to call the sine/cosine methods of the class,
# otherwise resort to math.cos/sin.
def cos(arg):
  try:
    return arg.cos()
  except TypeError:
    return _math.cos(arg)

def sin(arg):
  try:
    return arg.sin()
  except TypeError:
    return _math.sin(arg)

# Lift copy function to top level namespace.
def copy(arg):
  return _copy.copy(arg)
