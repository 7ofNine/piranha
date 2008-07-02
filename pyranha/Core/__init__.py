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

try:
  import numpy as __numpy
  import matplotlib as __matplotlib
except:
  print "Warning: many of Pyranha's capabilities rely on numpy and matplotlib. Please consider installing these packages."

def copy(arg):
	"""Standard copy function. Lifted from the copy module."""
	import copy as __copy
	return __copy.copy(arg)

psym_manager = _Core.__psym_manager()
expo_truncator = _Core.__expo_truncator()
norm_truncator = _Core.__norm_truncator()

def gui():
  try:
    import pyranha.Gui
    pyranha.Gui.mw.show()
  except ImportError:
    print "Gui support was not built."

def tc(args, flambda, t0 = 0., t1 = None, step = None):
  try:
    args_list = list(args)
  except TypeError:
    args_list = [args]
  if (t1 == None and step != None) or (t1 != None and step == None):
    raise TypeError, 'You must specify both a final time AND a step size.'
  if t1 == None:
    args_list_eval = map(lambda x: x.eval(t0),args_list)
    return abs(flambda(*args_list).eval(t0)-flambda(*args_list_eval))
  else:
    if (t1 <= t0):
      raise ValueError, 't1 must be strictly greater than t0.'
    time_array = __numpy.arange(t0,t1,step)
    retval = __numpy.array([],dtype=float)
    retval.resize(len(time_array))
    result = flambda(*args_list);
    j=0
    for t in time_array:
      args_list_eval = map(lambda x: x.eval(t),args_list)
      retval[j]=abs(result.eval(t)-flambda(*args_list_eval))
      j=j+1
    return time_array,retval;

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
							self.__args_descr = series.__arguments_description__.split()
							self.__args_list = series.__arguments__.__repr__().split(\"\\n\")
							self.__args_list.pop()
						def __get_tuple__(self):
							tmp_arg_list = [[int(x.split(None,1)[0]),x.split(None,1)[1]] for x in self.__args_list]
							retval = []
							n = 0
							prev = None
							for i in tmp_arg_list:
								if i[0] <= prev:
									n+=1
								tmp = [self.__args_descr[n]]
								tmp += i
								retval.append(tuple(tmp))
								prev = i[0]
							return tuple(retval)
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
