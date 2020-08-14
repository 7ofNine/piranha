# -*- coding: iso-8859-1 -*-
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

def check_order(s1,s2,n):
    l = len(s1.arguments)
    assert(l == len(s2.arguments))
    print(" check_order: " + str(n) + "\\n")
    return s1.filtered([lambda t: (t[0] * t[1]).order() < n] * l) == s2.filtered([lambda t: (t[0] * t[1]).order() < n] * l)
