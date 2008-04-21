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

# Handy definitions of common mathematical functions: try to call the sine/cosine methods of the class,
# otherwise resort to math.cos/sin.
def cos(arg):
  try:
    return arg.cos()
  except TypeError:
    return __math.cos(arg)

def sin(arg):
  try:
    return arg.sin()
  except TypeError:
    return __math.sin(arg)

# Lift copy function to top level namespace.
def copy(arg):
  return __copy.copy(arg)

psym_manager = _Core.__psym_manager()
expo_truncator = _Core.__expo_truncator()
norm_truncator = _Core.__norm_truncator()
