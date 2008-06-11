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
__type_list = []
for __i in __manipulators__:
	exec "import %s as __cur_manip" % __i
	__type_list.append(getattr(__cur_manip,__i.lower()))
	# Let's try to see if we can get the complex counterpart.
	try:
		__type_list.append(getattr(__cur_manip,__i.lower()+'c'))
	except:
		pass

manipulators_type_tuple = tuple(__type_list)

# Let's add the method to retrieve an index tuple to the series.
def __series_get_indices_tuple(self):
	retval = []
	i = 0
	go = True
	while go:
		try:
			tmp_index = getattr(self,'index'+str(i))
			retval.append(tmp_index)
			i += 1
		except AttributeError:
			go = False
	return tuple(retval)

for __i in __manipulators__:
	exec "import %s as __cur_manip" % __i
	exec("__cur_manip.%s.indices_tuple = __series_get_indices_tuple" % __i.lower())
	# Let's try to see if we can get the complex counterpart.
	try:
		exec("__cur_manip.%s.indices_tuple = __series_get_indices_tuple" % (__i.lower()+'c'))
	except:
		pass

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
