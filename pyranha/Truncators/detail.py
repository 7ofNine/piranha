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

class __truncators(object):
	def __init__(self):
		from pyranha.Truncators import _Truncators
		self.__list = filter(lambda x: x.endswith('_truncator'),dir(_Truncators))
		for n in self.__list:
			exec('self.%s = _Truncators.%s()' % (n.split('_truncator')[0][2:],n))
	def __repr__(self):
		from pyranha.Truncators import _Truncators
		retval = ''
		for n in self.__list:
			exec('t = _Truncators.%s()' % n)
			retval += ('%s: %s\n' % (n.split('_truncator')[0][2:],str(t)))
		return retval
	def unset(self):
		from pyranha.Truncators import _Truncators
		_Truncators.unset()

