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

class __psyms(object):
	def __repr__(self):
		from pyranha.Core import Psym
		l = Psym.list()
		retval = ''
		if not l:
			return 'No symbols defined.'
		max_length = max([len(s.name) for s in l])
		return reduce(lambda a,b: a + b, ['Symbol: \'' + s.name + '\'' + (' ' * (max_length - len(s.name))) + ' - ' + str(list(s.time_eval)) + '\n' for s in l])[0:-1]
	def __call__(self,names):
		"""
		Create psyms from a string of space-separated names.
		"""
		try:
			import IPython.ipapi
		except ImportError:
			raise ImportError("IPython not available.")
		ip = IPython.ipapi.get()
		for i in names.split():
			try:
				# Try to fetch the psym from the psym manager.
				ip.ex("%s = psym(\"%s\")" % (i,i))
			except SyntaxError:
				raise SyntaxError("The name '" + i + "' is not valid Python syntax, skipping.")

class __series(object):
	def __repr__(self):
		try:
			import IPython.ipapi
		except ImportError:
			raise ImportError("IPython not available.")
		import pyranha
		ip_ns = IPython.ipapi.get().user_ns
		l = [(
				x,
				ip_ns[x].__short_type__,
				len(ip_ns[x])
			) for x in [x for x in ip_ns if type(ip_ns[x]) in pyranha.manipulators]]
		l.sort()
		if not l:
			return 'No series defined.'
		m0, m1 = max([len(t[0]) for t in l]), max([len(t[1]) for t in l])
		return reduce(lambda a,b: a+ b,
			['Series: \'' + t[0] + '\'' + (' ' * (m0 - len(t[0]))) + ' - \'' + t[1] + '\'' + (' ' * (m1 - len(t[1]))) + ' - ' + str(t[2]) + '\n' for t in l])[0:-1]
	def __call__(self,names,series_t = None):
		"""
		Create series from a string of space-separated names. If the optional parameter series_t is None,
		Pyranha's default series type (ds) will be used. Otherwise the type series_t is used for
		series creation.
		"""
		try:
			import IPython.ipapi
		except ImportError:
			raise ImportError("IPython not available.")
		import pyranha
		if series_t == None:
			s_type = "ds"
		else:
			if series_t not in pyranha.manipulators:
				raise TypeError(str(type(series_t)) + " is not recognized as a valid series type.")
			s_type = series_t().__short_type__
		ip = IPython.ipapi.get()
		for i in names.split():
			try:
				# Try to fetch the psym from the psym manager.
				ip.ex("%s = %s(psym(\"%s\"))" % (i,s_type,i))
			except SyntaxError:
				raise SyntaxError("The name '" + i + "' is not valid Python syntax, skipping.")
