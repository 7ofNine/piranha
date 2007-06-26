# Copyright (C) 2007 by Francesco Biscani
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
from copy import *

# Global variables: theories of motion

#global elp1
#global elp2
#global elp3
#elp1 = np("elp1.csv")
#elp2 = (math.pi/2.) - np("elp2.csv")
#elp3 = np("elp3.csv")

def plot_ps(ps,mark='o',color='b',logscale=True):
	rng=range(ps.length())
	amps=list()
	for i in ps:
		amps.append(i.norm())
        semilogy(rng,amps,color+mark)
        xlim(0,len(ps))


def plot_ts(ps,t0,t1,n,xlab="Time",ylab="Value",mark='o',color='b'):
        rng=range(n)
        y=list()
        x=list()
        t=t0;
        delta=(t1-t0-0.)/n
        for i in rng:
            x.append(t)
            y.append(ps.t_eval(t))
            t+=delta
        plot(x,y,color+mark)
        xlabel(xlab)
        ylabel(ylab)

def plot_tc(tc,xlab="Time",ylab="Value",mark='+',logscale=True,lgnd=True):
	if logscale:
		plot_func=semilogy
	else:
		plot_func=plot
	rng=range(tc.size())
	time=list()
	exact=list()
	diff=list()
	for i in rng:
		time.append(tc.time(i))
		exact.append(abs(tc.hs(i)))
		diff.append(tc.error(i))
	plot_func(time,diff,'g'+mark,time,exact,'r',linewidth=2)
	if lgnd:
		legend(('Diff','Exact'),loc='best')
	xlabel(xlab)
	ylabel(ylab)
	title(r'$ \sigma = '+str(tc.sigma())+' $')

def plot_sc(sc,xlab="Term Index",ylab="Coefficient Delta"):
	rng=range(sc.size())
	col=list()
	for i in rng:
		col.append(sc.diffs(i))
	width=1.
	xlocations=array(range(len(col)))+0.5
	bar(xlocations,col,width=width)
	xlabel(xlab)
	if sc.is_relative():
		ylabel("Relative "+ylab)
	else:
		ylabel(ylab)
	xlim(0, xlocations[-1]+2*width)

def deg2rad(degrees,minutes,seconds):
	conv_ratio = math.pi/180.
	return (degrees*conv_ratio+minutes/60.*conv_ratio+seconds/3600.*conv_ratio)

def elp2000(time):
	return [sph_to_x(elp3.t_eval(time),elp2.t_eval(time),elp1.t_eval(time)),
		sph_to_y(elp3.t_eval(time),elp2.t_eval(time),elp1.t_eval(time)),
		sph_to_z(elp3.t_eval(time),elp2.t_eval(time),elp1.t_eval(time))]

def doodson_Bnm(n,m):
	if n < 2:
		print "Invalid n: it must be >= 2"
		return npc()
	return (wig_rot(n,m,0,-astro.eps_0(),0,elp2,elp1)*elp3.complex_multiangle(0,-m))*natural_pow(n+1,elp3.inv())

def benchmark(filename="prec_test",steps=1000):
	foo = wig_rot(1,1,1.,2.,3.,elp2,elp1)
	print "Final length=",foo.real().length(),",",foo.imag().length()
	bench=tc_wig_rot_np(foo,0.,.01,steps,1,1,1.,2.,3.,elp2,elp1)
	bench.gnuplot_save(filename);
	print "Pack RATIO: ",pack_ratio()
	return bench

def tidal_accel(m,d,r):
	a_A=abs(astro.G()*m*((2.*d*r-r*r)/(d*d*(d-r)*(d-r))))
	a_B=abs(astro.G()*m*((-2.*d*r-r*r)/(d*d*(d+r)*(d+r))))
	return a_A, a_B
