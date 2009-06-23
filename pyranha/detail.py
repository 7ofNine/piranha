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

def __series_filtered(self, criterion = None):
	"""
	Return input series with terms filtered according to criterion.

	If criterion is a callable binary object, the series will be decomposed into a list of coefficient-key pairs and the callable
	will be invoked using coefficients and keys as arguments. If the return value of the callable is True, then
	coefficient and key are retained, otherwise they are eliminated from the series.

	If criterion is a list of callable binary objects, each callable will be called recursively to filter out
	coefficient-key pairs, beginning with the current series and descending in a recursive fashion into the coefficient
	series (and from there into the coefficient series of the coefficient series, etc.).

	A None criterion is interpreted as a callable that always returns True.

	Please note that since reassembling coefficient-key pairs involves series multiplications, active truncators
	have an effect on the filtering process.
	"""
	from copy import copy
	if criterion is None:
		return copy(self)
	try:
		iter(criterion)
		crit = list(criterion)
	except TypeError:
		crit = [criterion]
	# If list is empty, return self.
	if not crit:
		return copy(self)
	retval = type(self)()
	rec_depth = len(crit)
	if rec_depth > len(self.arguments):
		raise ValueError('Cannot apply %d recursive filters to a series of echelon level %d' % (rec_depth,len(self.arguments) - 1))
	for c in crit:
		if not callable(c) and not c is None:
			raise ValueError('Please provide binary callables (or None) as filtering criterions.')
	def filter_series(cur,tot,s,crits):
		if cur > tot:
			return copy(s)
		retval = type(s)()
		l = s.split(cur)
		for t in l:
			if crits[cur] is None or crits[cur](t[0],t[1]):
				retval += filter_series(cur + 1, tot, t[0], crits) * t[1]
		return retval
	return filter_series(0,rec_depth - 1,self,crit)

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

def __series_split(self, n = 0):
	"""
	Split the series into a sequence of coefficient-key pairs.

	For n > 0, this method will try to split the series at higher echelon levels. That is,
	if n == 1 and the series is degenerate (i.e., one single term with unitary key), a split()
	on the single term's coefficient series will take place and the return value will be a
	sequence of coefficient-key pairs for the coefficient series.
	"""
	return self.__split__(n)

def __series_eval_sub(self,subs):
	"""
	Evaluate by substitution.
	"""
	from copy import copy
	from pyranha.Core import psym
	args_set = list(set(reduce(lambda x,y: list(x) + list(y), self.arguments)))
	retval = copy(self)
	for a in args_set:
		if not a.name in subs:
			raise ValueError('The provided substitution dictionary does not contain all the arguments of the series.')
		retval = retval.sub(psym(a.name),type(self)(subs[a.name]))
	return retval.eval(0)

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
		__add_method(i, "split", __series_split)
		__add_method(i, "eval_sub", __series_eval_sub)

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
