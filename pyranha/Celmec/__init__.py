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

import pyranha as __pyranha

def r_a(e,M):
  """
  Calculate the elliptic expansion of r/a in terms of e and M.
  """
  try:
    return __pyranha.ds.r_a(e,M)
  except AttributeError:
    raise AttributeError, "The series type '" + str(ds) +  "' does not offer an r_a method."

def sin_f(e,M):
  """
  Calculate the elliptic expansion of sin(f) in terms of e and M.
  """
  try:
    return __pyranha.ds.sin_f(e,M)
  except AttributeError:
    raise AttributeError, "The series type '" + str(ds) +  "' does not offer a sin_f method."

def cos_E(e,M):
  """
  Calculate the elliptic expansion of cos(E) in terms of e and M.
  """
  try:
    return __pyranha.ds.cos_E(e,M)
  except AttributeError:
    raise AttributeError, "The series type '" + str(ds) +  "' does not offer a cos_E method."

def E(e,M):
  """
  Calculate the elliptic expansion of E in terms of e and M.
  """
  try:
    return __pyranha.ds.E(e,M)
  except AttributeError:
    raise AttributeError, "The series type '" + str(ds) +  "' does not offer an E method."
