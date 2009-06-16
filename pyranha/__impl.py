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

def __series_filtered(series, criterion = None):
	"""
	Return input series with terms filtered according to criterion.
	"""
	if criterion is None:
		return series
	if not callable(criterion):
		raise ValueError("You must provide a unary callable object for series filtering.")
	retval = type(series)()
	for t in series:
		if criterion(t):
			retval += t[0] * t[1]
	return retval

def __build_manipulators(manipulators):
	"""
	Build the list of manipulator types
	"""
	type_list = []
	for i in manipulators:
		exec "import pyranha.%s as cur_manip" % i
		type_list.append(getattr(cur_manip,i.lower()))
		# Let's try to see if we can get the complex counterpart.
		try:
			type_list.append(getattr(cur_manip,i.lower()+'c'))
		except AttributeError:
			pass
	return tuple(type_list)

def __series_short_type(self):
	"""
	Return a short string containing the series' class name.
	"""
	return str(type(self)).rpartition('.')[-1].strip('>\'')

def __series_psi(self, n = 0, s = 1):
	"""
	Return the limit of the power series expansion compatible with the truncator currently in use by the series.
	The optional arguments are the starting degree (n) and the step size (s) of the power series expansion.
	"""
	return self.__psi__(n,s)

def __add_method(module_name,method_name,function):
	"""
	Add a method to a manipulator.
	"""
	exec "import %s as __cur_manip" % module_name
	exec("__cur_manip.%s.%s = function" % (module_name.lower(),method_name))
	# Try to take care of the complex counterpart.
	try:
		exec("__cur_manip.%s.%s = function " % ((module_name.lower()+'c'),method_name))
	except AttributeError:
		pass

def __add_property(module_name, property_name, fget=None, fset=None, fdel=None, doc=None):
	"""
	Add a property to a manipulator.
	"""
	exec "import %s as __cur_manip" % module_name
	exec("__cur_manip.%s.%s = property(fget,fset,fdel,doc)" % (module_name.lower(),property_name))
	try:
		exec("__cur_manip.%s.%s = property(fget,fset,fdel,doc) " % ((module_name.lower()+'c'),property_name))
	except AttributeError:
		pass

def __enhance_manipulators(manipulators):
	"""
	Add useful methods and properties to a list of manipulators.
	"""
	for i in manipulators:
		__add_property(i, "__short_type__", __series_short_type)
		__add_method(i, "filtered", __series_filtered)
		__add_method(i, "psi", __series_psi)

import pyranha

# Let's try do define a default series type
if len(filter(lambda x: x not in pyranha.__manipulators__,pyranha.__all__)) > 0:
	exec "import %s as __last_manipulator" % pyranha.__manipulators__[-1]
	ds = getattr(__last_manipulator,pyranha.__manipulators__[-1].lower());
	setattr(pyranha,'ds',ds)
	print "Default series type is " + str(ds)
else:
	print "Default series type could not be established, assigning None."
	ds = None

__enhance_manipulators(pyranha.__manipulators__)

manipulators = tuple(__build_manipulators(pyranha.__manipulators__))
