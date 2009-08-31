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

__doc__ = 'Keplerian motion perturbed by constant thrust.'

axis = 0
thrust = 0.01

def H1(mdelaunay):
	"""
	Perturbed hamiltonian.
	"""
	from pyranha.Celmec import mdelaunay2oe, cos_E, sin_E, orbitalR
	from pyranha.Math import root, cos, sin
	from numpy import dot
	# TODO: check input arguments.
	# Calculate classical orbital elements from input mod delaunay elements.
	a, e, s, omega, Omega, M = mdelaunay2oe(mdelaunay)
	# Position vector in the orbital plane.
	ro = [a * (cos_E(e,M) - e),a * root(2,1 - e ** 2) * sin_E(e,M), 0]
	# Build the orbital rotation matrix.
	Rxq = orbitalR([cos(omega),sin(omega),1 - 2 * s **2,2 * s * root(2,1 - s ** 2),cos(Omega),sin(Omega)])
	# Position vector in space.
	r = dot(Rxq,ro)
	return -r[axis]

def numerical(s0,tf,n):
	"""
	Return a list of n state vectors starting from t = 0 to t = tf assuming as initial conditions
	the state s0 (list of 6 numbers, rectangular position + velocity).
	"""
	def dyn(s,_):
		from numpy.linalg import norm
		r = norm(s[:3])
		r3 = r * r * r
		retval = [s[3],s[4],s[5],-s[0] / r3,-s[1] / r3,-s[2] / r3]
		retval[3 + axis] = retval[3 + axis] + thrust
		return retval
	from scipy.integrate import odeint
	from numpy import linspace
	return odeint(dyn,s0,linspace(0,tf,n))

def s2oe(s0):
	from numpy.linalg import norm
	from numpy import array, cross, dot
	from math import acos, atan2, cos, pi, sin, sqrt
	# Position and velocity.
	v_r = s0[:3]
	v_v = s0[3:]
	r = norm(v_r)
	v = norm(v_v)
	# Angular momentum.
	v_h = cross(v_r,v_v)
	h = norm(v_h)
	# Semi-major axis.
	a = -1. / (2. * (v * v / 2. - 1. / r))
	# Eccentricity vector.
	v_e = cross(v_v,v_h) - array(v_r) / r
	e = norm(v_e)
	# s = sin(i/2)
	s = sqrt((1. - v_h[2] / h) / 2.)
	# Vector of the ascending node.
	v_n = [-v_h[1],v_h[0],0]
	n = norm(v_n)
	# omega and Omega.
	if e == 0:
		om = 0
	elif s == 0:
		om = acos(v_e[0] / e)
	elif v_e[2] >= 0:
		om = acos(dot(v_n,v_e) / (n * e))
	else:
		om = 2 * pi - acos(dot(v_n,v_e) / (n * e))
	if s == 0:
		Om = 0
	elif v_n[1] >= 0:
		Om = acos(v_n[0] / n)
	else:
		Om = 2 * pi - acos(v_n[0] / n)
	# True anomaly.
	if e == 0 and s == 0:
		if v_v[0] <= 0:
			f = acos(v_r[0] / r)
		else:
			f = 2 * pi - acos(v_r[0] / r)
	elif e == 0:
		if dot(v_n,v_v) <= 0:
			f = acos(dot(v_n,v_r) / (n * r))
		else:
			f = 2 * pi - acos(dot(v_n,v_r) / (n * r))
	elif dot(v_r,v_v) >= 0:
		f = acos(dot(v_e,v_r) / (e * r))
	else:
		f = 2 * pi - acos(dot(v_e,v_r) / (e * r))
	# Eccentric anomaly.
	E = atan2(sqrt(1 - e * e) * sin(f),e + cos(f))
	# Mean anomaly.
	M = E - e * sin(E)
	return [a,e,s,om,Om,M]
