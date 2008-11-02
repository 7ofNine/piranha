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
	"""Standard copy function. Lifted from the copy module."""
	import copy as __copy
	return __copy.copy(arg)

def __create_psyms(self, names):
	"""
	Create psyms from a string of space-separated names. If a name is already assigned to a psym,
	create a copy of the psym with the same properties of its homonymous.
	"""
	try:
		import IPython.ipapi
	except ImportError:
		raise ImportError("IPython not available.")
	ip = IPython.ipapi.get()
	for i in names.split():
		try:
			# Try to fetch the psym from the psym manager.
			ip.ex("%s = psyms[\"%s\"]" % (i,i))
		except UserWarning:
			# If we failed, create a new psym.
			ip.ex("%s = psym(\"%s\")" % (i,i))
		except SyntaxError:
			print "The name '" + i + "' is not valid Python syntax, skipping."

setattr(_Core.__psyms,"__call__",__create_psyms)
psyms = _Core.__psyms()

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
			ip.ex("%s = %s(psyms[\"%s\"])" % (i,s_type,i))
		except UserWarning:
			# If we failed, create a new psym.
			ip.ex("%s = %s(psym(\"%s\"))" % (i,s_type,i))
		except SyntaxError:
			print "The name '" + i + "' is not valid Python syntax, skipping."

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

degree_truncator = _Core.__degree_truncator()
expo_truncator = _Core.__expo_truncator()
norm_truncator = _Core.__norm_truncator()
settings = _Core.__settings()

def gui():
	try:
		import pyranha.Gui
		pyranha.Gui.mw.show()
	except ImportError:
		print "Gui support is not available or PyQt4 is not installed."

class tc(object):
	def __init__(self, args, f, t0, t1, step, res = None):
		try:
			import numpy
		except ImportError:
			raise ImportError("Numpy is not available.")
		# If args is a tuple, transform it into a list, if args is a single value
		# build a list from it.
		if isinstance(args,type(())):
			args_list = list(args)
		else:
			args_list = [args]
		if t1 <= t0:
			raise ValueError, 't1 must be strictly greater than t0.'
		if step <= 0:
			raise ValueError, 'Step must be strictly positive.'
		if res != None:
			r = res
		else:
			r = f(*args_list)
		args_list_eval_funcs = map(lambda x: x.eval,args_list)
		self.time_array = numpy.arange(t0,t1,step)
		self.eval_array = numpy.array(map(lambda t: r.eval(t), self.time_array))
		self.diff_array = numpy.array(map(lambda x: abs(x[1] - f(*map(lambda y: y(x[0]),args_list_eval_funcs))),zip(self.time_array,self.eval_array)))
	def plot(self):
		import pylab
		pylab.semilogy(self.time_array,self.eval_array,self.time_array,self.diff_array)

def __decorate_args_tuples():
	import _Core
	n = 1
	try:
		while True:
			base_class_name = "_Core.__base_args_tuple"+str(n)+"__"
			new_class_name = "__args_tuple"+str(n)+"__"
			repr_function_name = "__args_tuple"+str(n)+"_repr__"
			exec(
				"""class %s(%s):
						def __init__(self,series):
							%s.__init__(self,series.__arguments__)
							self.__args_descr = series.__arguments_description__.split(\"\\n\")
							self.__args_list = series.__arguments__.__repr__().split(\"\\n\")
							self.__args_descr.pop()
							self.__args_list.pop()
						def __get_tuple__(self):
							tmp = []
							for i in self.__args_descr:
								tmp.extend([i.split(":")[0]]*int(i.split(":")[1]))
							return tuple([(i[0],i[1].split()[0],i[1].split()[1]) for i in zip(tmp,self.__args_list)])
						def __repr__(self):
							retval = \"\"
							for i in self.__get_tuple__():
								retval += str(i[0])+\" \"+str(i[1])+\" \"+str(i[2])+\"\\n\"
							return retval"""
				% (new_class_name, base_class_name, base_class_name))
			exec("setattr(_Core, new_class_name, %s)" % new_class_name)
			n += 1
	except AttributeError:
		pass

__decorate_args_tuples()
