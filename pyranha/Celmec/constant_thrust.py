# -*- coding: utf-8 -*-
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
thrust = 0.001

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

def build_analytical(s0,pert_order):
	from copy import copy
	from pyranha import truncators
	from pyranha.Celmec import oe2mdelaunay, lieS
	from pyranha.Math import integrate, partial
	from pyranha.Core import psym, rational
	from pyranha.Qqps import qqps
	import pyranha.Celmec.constant_thrust as ct
	# Some parameters.
	# NOTE: these could depend on a parametric numerical precision.
	H1_exp_order = 2 # Order of ecc + incl in the initial expansion of the perturbing Hamiltonian.
	mom = ['Lam','P','Q']
	coord = ['lam','p','q']
	# Create needed series.
	Lam, P, Q, lam, p, q, lam0, eps = [qqps(psym(n)) for n in mom + coord + ['lam0','eps']]
	# Set truncators for initial expansion of perturbing Hamiltonian.
	truncators.unset()
	truncators.degree.set(['P','Q'],rational(H1_exp_order + 1,2))
	# Hamiltonian.
	H = -1 * (2 * Lam ** 2) ** -1 + eps * H1([Lam,P,Q,lam,p,q])
	# Solver for the homological equation.
	def he_s(Hk):
		# This will integrate wrt lambda the part of the k-th order of the Hamiltonian that depends on the coordinates
		# and return the result as the generator for the Lie transform.
		# return integrate(Lam ** 3 * Hk.filtered([lambda t: any(n in t[0] or n in t[1] for n in coord)]),'lam')
		return integrate(Lam ** 3 * Hk.filtered(lambda t: t[1].h_order('lam') != 0),'lam')
	# Clear return values.
	# Replace lambda with delta_lam + lambda_0.
	ct.H_list = [H.sub('Q',qqps())]
	ct.chi_list = []
	ct.direct = [[Lam, P, Q, lam, p, q]]
	ct.inverse = [[Lam, P, Q, lam, p, q]]
	ct.freqs = []
	# Start the iteration.
	for i in range(0,pert_order):
		cur_H = ct.H_list[-1]
		print('\033[1;32mCurrent perturbative order: %d\033[1;m' % (i + 1))
		Hk = cur_H.filtered([None,lambda t: t[1].degree('eps') == i + 1]) * eps ** (-i -1) 
		print('Solving homological equation...')
		chi = he_s(Hk)
		print('Calculating Hamiltonian...')
		# Truncators setup.
		truncators.degree.set(['eps'],pert_order + 1)
		l = eps.psi()
		truncators.degree.set(['P','Q'],rational(H1_exp_order + 1,2))
		ct.H_list.append(lieS(eps ** (i + 1),chi,cur_H,mom,coord,l))
		ct.chi_list.append(chi)
		print('Calculating ' +'\033[1;31mdirect\033[1;m ' +'coordinate transformation...')
		ct.direct.append([lieS(eps ** (i + 1),chi,arg,mom,coord,l) for arg in ct.direct[-1]])
		print('Calculating ' +'\033[1;31minverse\033[1;m ' +'coordinate transformation...')
		ct.inverse.append([lieS(eps ** (i + 1),-chi,arg,mom,coord,l) for arg in ct.inverse[-1]])
	# Build dictionaries for evaluation.
	var_names = mom + coord + ['eps','two','lam0']
	eval_direct = dict(zip(var_names,oe2mdelaunay(s2oe(s0)) + [ct.thrust,2.,0]))
	# We have to take into account the modification of lambda for delta_lam.
	eval_direct['lam0'] = eval_direct['lam']
	eval_direct['lam'] = 0
	# Now te initial conditions for the inverse.
	eval_inverse = dict(zip(var_names,[ct.inverse[-1][i].eval(eval_direct) for i in range(0,6)] + [ct.thrust,2.,eval_direct['lam0']]))
	# Now the frequencies of the angles.
	ct.freqs = dict(zip(coord,[partial(ct.H_list[-1],arg).eval(eval_inverse) for arg in mom]))
	def t_eval(t):
		retval = copy(eval_inverse)
		for n in coord:
			retval[n] += freqs[n] * t
		return retval
	ct.t_eval = t_eval

def build_poincare(s0):
	from copy import copy
	from pyranha import truncators
	from pyranha.Celmec import oe2mdelaunay, mdelaunay2poincare, lieS
	from pyranha.Math import integrate, partial
	from pyranha.Core import psym, rational
	from pyranha.Qqps import qqps, qqpsc
	import pyranha.Celmec.constant_thrust as ct
	# Some parameters.
	# NOTE: these could depend on a parametric numerical precision.
	pert_order = 1
	H1_exp_order = 4 # Order of ecc + incl in the initial expansion of the perturbing Hamiltonian.
	md_mom = ['Lam','P','Q']
	md_coord = ['lam','p','q']
	# Create needed series.
	Lam, P, Q, lam, p, q, lam0, eps = [qqps(psym(n)) for n in md_mom + md_coord + ['lam0','eps']]
	# Set truncators for initial expansion of perturbing Hamiltonian.
	truncators.unset()
	truncators.degree.set(['P','Q'],rational(H1_exp_order + 1,2))
	# Hamiltonian.
	H = -1 * (2 * Lam ** 2) ** -1 + eps * H1([Lam,P,Q,lam,p,q])
	truncators.unset()
	# Transform into Poincare' variables.
	p_mom = ['Lam','y','z']
	p_coord = ['lam','x','v']
	Lam, y, z, lam, x, v, two = [qqps(psym(n)) for n in p_mom + p_coord + ['two']]
	H = H.ei_sub('p',qqpsc(y * (two * P) ** rational(-1,2),x * (two * P) ** rational(-1,2))).ei_sub('q',qqpsc(z * (two * Q) ** rational(-1,2),v * (two * Q) ** rational(-1,2))).sub('two',qqps(2)).sub('P',(x ** 2 + y ** 2) / 2).sub('Q',(z ** 2 + v ** 2) / 2)
	# Solver for the homological equation.
	def he_s(Hk):
		# This will integrate wrt lambda the part of the k-th order of the Hamiltonian that depends on the coordinates
		# and return the result as the generator for the Lie transform.
		return integrate(Lam ** 3 * Hk.filtered([lambda t: any(n in t[0] or n in t[1] for n in p_coord)]),'lam')
		#return integrate(Lam ** 3 * Hk.filtered(lambda t: t[1].h_order('lam') != 0),'lam')
	# Truncator for Lie Series expansions.
	truncators.degree.set(['eps'],pert_order + 1)
	# Clear return values.
	# Replace lambda with delta_lam + lambda_0.
	ct.H_list = [H.sub('lam',lam + lam0)]
	ct.chi_list = []
	ct.direct = [[Lam, y, z, lam, x, v]]
	ct.inverse = [[Lam, y, z, lam, x, v]]
	ct.freqs = []
	# Start the iteration.
	for i in range(0,pert_order):
		cur_H = ct.H_list[-1]
		print('\033[1;32mCurrent perturbative order: %d\033[1;m' % (i + 1))
		Hk = cur_H.filtered([None,lambda t: t[1].degree('eps') == i + 1]) * eps ** (-i -1) 
		print('Solving homological equation...')
		chi = he_s(Hk)
		print('Calculating Hamiltonian...')
		ct.H_list.append(lieS(eps ** (i + 1),chi,cur_H,p_mom,p_coord))
		ct.chi_list.append(chi)
		print('Calculating ' +'\033[1;31mdirect\033[1;m ' +'coordinate transformation...')
		ct.direct.append([lieS(eps ** (i + 1),chi,arg,p_mom,p_coord) for arg in ct.direct[-1]])
		print('Calculating ' +'\033[1;31minverse\033[1;m ' +'coordinate transformation...')
		ct.inverse.append([lieS(eps ** (i + 1),-chi,arg,p_mom,p_coord) for arg in ct.inverse[-1]])
	# Build dictionaries for evaluation.
	var_names = p_mom + p_coord + ['eps','two','lam0']
	eval_direct = dict(zip(var_names,mdelaunay2poincare(oe2mdelaunay(s2oe(s0))) + [ct.thrust,2.,0]))
	# We have to take into account the modification of lambda for delta_lam.
	eval_direct['lam0'] = eval_direct['lam']
	eval_direct['lam'] = 0
	# Now te initial conditions for the inverse.
	eval_inverse = dict(zip(var_names,[ct.inverse[-1][i].eval(eval_direct) for i in range(0,6)] + [ct.thrust,2.,eval_direct['lam0']]))
	# Now the frequencies of the angles.
	ct.freqs = dict(zip(p_coord,[partial(ct.H_list[-1],arg).eval(eval_inverse) for arg in p_mom]))
	def t_eval(t):
		retval = copy(eval_inverse)
		for n in p_coord:
			retval[n] += freqs[n] * t
		return retval
	ct.t_eval = t_eval

def s2oe(s):
	from numpy.linalg import norm
	from numpy import array, cross, dot
	from math import acos, atan2, cos, pi, sin, sqrt
	# Position and velocity.
	v_r = s[:3]
	v_v = s[3:]
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
		om = 0.
	elif s == 0:
		om = acos(v_e[0] / e)
	elif v_e[2] >= 0:
		om = acos(dot(v_n,v_e) / (n * e))
	else:
		om = 2 * pi - acos(dot(v_n,v_e) / (n * e))
	if s == 0:
		Om = 0.
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
