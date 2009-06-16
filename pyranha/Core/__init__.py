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

from _Core import *

def copy(arg):
	"""
	Standard copy function. Lifted from the copy module.
	"""
	import copy as __copy
	return __copy.copy(arg)

def latex(arg):
	"""
	Return latex-formatted string representation of the object.
	"""
	if not hasattr(arg,"_latex_"):
		raise AttributeError('Object does not provide a _latex_() method.')
	return arg._latex_()

def psyms(names):
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

def series(names,series_t = None):
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
		if series_t not in pyranha.manipulators_type_tuple:
			raise TypeError(str(type(series_t)) + " is not recognized as a valid series type.")
		s_type = series_t().__short_type__
	ip = IPython.ipapi.get()
	for i in names.split():
		try:
			# Try to fetch the psym from the psym manager.
			ip.ex("%s = %s(psym(\"%s\"))" % (i,s_type,i))
		except SyntaxError:
			raise SyntaxError("The name '" + i + "' is not valid Python syntax, skipping.")

def load(*args):
	"""
	Load series from list of filenames arguments.
	Type of series will be inferred - if possible - from the files' extensions.
	"""
	try:
		import IPython.ipapi
	except ImportError:
		raise ImportError("IPython not available.")
	ip = IPython.ipapi.get()
	args_list = list(args)
	for i in args_list:
		s = i.split(".")
		var_name = s.pop(0)
		if not s:
			raise ValueError("Could not establish an extension for the filename \"" + arg + "\".")
			return
		ext_name = s.pop()
		try:
			ip.ex("%s = %s(\"%s\")" % (var_name,ext_name,i))
		except NameError:
			raise NameError("Series type \"" + ext_name + "\" is unknown.")

settings = _Core.__settings()

def gui():
	try:
		import pyranha.Gui
		pyranha.Gui.mw.show()
	except ImportError:
		print "Gui support is not available or PyQt4 is not installed."

norm_truncator = _Core.__norm_truncator()
degree_truncator = _Core.__degree_truncator()
