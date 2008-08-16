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

#from numpy import *
#from matplotlib.pylab import *

__all__ = [@MODULE_LIST@]
__manipulators__ = [@MANIPULATOR_LIST@]

print "Pyranha initializing..."
print "Available manipulators: ", __manipulators__
print "Other modules: ", filter(lambda x: x not in __manipulators__,__all__)

# Let's try do define a default series type
if len(filter(lambda x: x not in __manipulators__,__all__)) > 0:
	exec "import %s as __last_manipulator" % __manipulators__[-1]
	ds = getattr(__last_manipulator,__manipulators__[-1].lower());
	print "Default series type is " + str(ds)
else:
	print "Default series type could not be established, assigning None."
	ds = None

# Let's build the list of manipulator types
def __build_manipulators_type_tuple():
	__type_list = []
	for i in __manipulators__:
		exec "import %s as __cur_manip" % i
		__type_list.append(getattr(__cur_manip,i.lower()))
		# Let's try to see if we can get the complex counterpart.
		try:
			__type_list.append(getattr(__cur_manip,i.lower()+'c'))
		except AttributeError:
			pass
	return tuple(__type_list)

manipulators_type_tuple = tuple(__build_manipulators_type_tuple())

def __series_short_type(self):
	return str(type(self)).rpartition('.')[-1].strip('>\'')

def __series_filter(self, *args):
	self.__set_shared_arguments__()
	new_series = type(self)()
	new_series.__set_arguments__(self.__arguments__)
	map(lambda t: new_series.__append__(t),__base_series_filter(self,*args))
	new_series.__trim__()
	return new_series

def __base_series_filter(series,*args):
	import copy
	if not hasattr(series,"__iter__"):
		# Input is not iteratable, i.e. it is a numerical coefficient. Retval will be a copy of input.
		return copy.copy(series)
	args_list = list(args)
	func = None
	if args_list:
		func = args_list[0]
		if func != None and not callable(func):
			raise ValueError("Please provide functions or None as arguments for filter().")
	else:
		# List of functions is exhausted, nothing to filter. Just use a copy of series as return value.
		return copy.copy(series)
	args_list.pop(0)
	retval = type(series)()
	for i in series:
		j = type(i)()
		j.key = i.key
		j.cf = __base_series_filter(i.cf,*args_list)
		if func == None or func(j):
			retval.__append__(j)
	return retval

def __series_arguments(self):
	import pyranha.Core._Core as _Core
	import re
	args_number = re.match('.*__base_args_tuple(..?)__.*',str(type(self.__arguments__))).group(1)
	exec("retval = _Core.__args_tuple%s__(self)" % args_number)
	return retval

def __series_plot(self, *args, **kwargs):
	"""Plot series' terms."""
	try:
		import pylab
	except ImportError:
		raise ImportError("Matplotlib is not available.")
	# Set the shared arguments.
	self.__set_shared_arguments__()
	# Let's start parsing the kwargs and plucking the non matplotlib ones.
	mpl_kw = kwargs
	try:
		log = mpl_kw.pop("log")
	except KeyError:
		log = False
	try:
		cmp = mpl_kw.pop("cmp")
	except KeyError:
		cmp = None
	try:
		key = mpl_kw.pop("key")
	except KeyError:
		key = None
	try:
		reverse = mpl_kw.pop("reverse")
	except KeyError:
		reverse = False
	try:
		f = mpl_kw.pop("f")
	except KeyError:
		# If f has not been specified in the arguments' list, use as f the key function, if
		# it is other than None. Otherwise use as f the norm of the term.
		if key == None:
			f = lambda t: t.cf.norm * t.key.norm
		else:
			f = key
	plot_iteratable = None
	if log:
		plot_func = pylab.semilogy
	else:
		plot_func = pylab.plot
	if key or cmp:
		plot_iteratable = sorted(self,cmp,key,reverse)
	else:
		plot_iteratable = self
	plot_iteratable = map(f,plot_iteratable)
	plot_func(range(len(plot_iteratable)),plot_iteratable,*args,**kwargs)

def __add_method(module_name,method_name,function):
	exec "import %s as __cur_manip" % module_name
	exec("__cur_manip.%s.%s = function" % (module_name.lower(),method_name))
	# Try to take care of the complex counterpart.
	try:
		exec("__cur_manip.%s.%s = function " % ((module_name.lower()+'c'),method_name))
	except AttributeError:
		pass

def __add_property(module_name, property_name, fget=None, fset=None, fdel=None, doc=None):
	exec "import %s as __cur_manip" % module_name
	exec("__cur_manip.%s.%s = property(fget,fset,fdel,doc)" % (module_name.lower(),property_name))
	try:
		exec("__cur_manip.%s.%s = property(fget,fset,fdel,doc) " % ((module_name.lower()+'c'),property_name))
	except AttributeError:
		pass

def __enhance_manipulators():
	for i in __manipulators__:
		__add_property(i, "__short_type__", __series_short_type)
		__add_method(i, "filter", __series_filter)
		__add_property(i, "arguments", __series_arguments)
		__add_method(i, "plot", __series_plot)

__enhance_manipulators()

print "Pyranha is ready."

# Global variables: theories of motion

#global elp1
#global elp2
#global elp3
#elp1 = np("elp1.csv")
#elp2 = (math.pi/2.) - np("elp2.csv")
#elp3 = np("elp3.csv")

#def plot_ps(ps,mark='o',color='b',logscale=True):
	#rng=range(ps.length())
	#amps=list()
	#for i in ps:
		#amps.append(i.norm())
        #semilogy(rng,amps,color+mark)
        #xlim(0,len(ps))

#def plot_sc(sc,xlab="Term Index",ylab="Coefficient Delta"):
	#rng=range(sc.size())
	#col=list()
	#for i in rng:
		#col.append(sc.diffs(i))
	#width=1.
	#xlocations=array(range(len(col)))+0.5
	#bar(xlocations,col,width=width)
	#xlabel(xlab)
	#if sc.is_relative():
		#ylabel("Relative "+ylab)
	#else:
		#ylabel(ylab)
	#xlim(0, xlocations[-1]+2*width)

#def deg2rad(degrees,minutes,seconds):
	#conv_ratio = math.pi/180.
	#return (degrees*conv_ratio+minutes/60.*conv_ratio+seconds/3600.*conv_ratio)

#def elp2000(time):
	#return [sph_to_x(elp3.t_eval(time),elp2.t_eval(time),elp1.t_eval(time)),
		#sph_to_y(elp3.t_eval(time),elp2.t_eval(time),elp1.t_eval(time)),
		#sph_to_z(elp3.t_eval(time),elp2.t_eval(time),elp1.t_eval(time))]

#def doodson_Bnm(n,m):
	#if n < 2:
		#print "Invalid n: it must be >= 2"
		#return npc()
	#return (wig_rot(n,m,0,-astro.eps_0(),0,elp2,elp1)*elp3.complex_multiangle(0,-m))*natural_pow(n+1,elp3.inv())

#def benchmark(filename="prec_test",steps=1000):
	#foo = wig_rot(1,1,1.,2.,3.,elp2,elp1)
	#print "Final length=",foo.real().length(),",",foo.imag().length()
	#bench=tc_wig_rot_np(foo,0.,.01,steps,1,1,1.,2.,3.,elp2,elp1)
	#bench.gnuplot_save(filename);
	#print "Pack RATIO: ",pack_ratio()
	#return bench

#def tidal_accel(m,d,r):
	#a_A=abs(astro.G()*m*((2.*d*r-r*r)/(d*d*(d-r)*(d-r))))
	#a_B=abs(astro.G()*m*((-2.*d*r-r*r)/(d*d*(d+r)*(d+r))))
	#return a_A, a_B
