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

#__all__ = [@MODULE_LIST@]
#__manipulators__ = [@MANIPULATOR_LIST@]
#
#print("Pyranha initializing...")
#print("Available manipulators: ", __manipulators__)
#print("Other modules: ", filter(lambda x: x not in __manipulators__,__all__))
#
#from detail import manipulators
#from Truncators import truncators
#
#print("Pyranha is ready.")

import os;
os.add_dll_directory("@PIRANHA_INSTALL_PATH@") # needed since 3.8 to be able to load dependency dll's

__all__ = [@MODULE_LIST@]
__manipulators__ = [@MANIPULATOR_LIST@]

print("Pyranha initializing...")
print("Available manipulators: ", __manipulators__)
print("Other modules: ", [x for x in __all__ if x not in __manipulators__])

from .detail import manipulators
from .Truncators import truncators

print("Pyranha is ready.")
