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

import pyranha.Celmec

class lt_divergent(pyranha.Celmec.lie_theory):
	def __init__(self,order,t,thrust,axis = 0):
		from pyranha.Core import psym
		from pyranha import manipulators
		if not t in manipulators:
			raise TypeError('t must be a series type.')
		if not isinstance(order,int) or order <= 0:
			raise ValueError('Invalid order.')
		if not isinstance(thrust,float) or thrust <= 0:
			raise ValueError('Invalid thrust.')
		if not isinstance(axis,int) or not axis in [0,1,2]:
			raise ValueError('Invalid axis.')
		self.__thrust = thrust
		self.__axis = axis
		Lam, P, Q, lam, p, q = t(psym('Lam')), t(psym('P')), t(psym('Q')), t(psym('lam')), t(psym('p')), t(psym('q'))
		lam0, dlam = t(psym('lam0')), t(psym('dlam'))
		eps_series = t(psym('eps'))
		H = - (2 * Lam ** 2) ** -1 + eps_series * self.__H1([Lam,P,Q,lam,p,q]).sub('lam',lam0 + dlam)
		super(lt_divergent,self).__init__(H,'eps',['Lam','P','Q'],['dlam','p','q'],[lambda Hn: Lam ** 3 * Hn.integrate('dlam')] * order)
	def __H1(self,mdelaunay):
		# Perturbed Hamiltonian.
		from pyranha.Core import rational
		from pyranha.Celmec import mdelaunay2oe, cos_E, sin_E, orbitalR
		from numpy import dot
		# Calculate classical orbital elements from input mod Delaunay elements.
		a, e, s, omega, Omega, M = mdelaunay2oe(mdelaunay)
		# Position vector in the orbital plane.
		ro = [a * (cos_E(e,M) - e),a * (1 - e ** 2) ** rational(1,2) * sin_E(e,M), 0]
		# Build the orbital rotation matrix.
		Rxq = orbitalR([omega.cos(),omega.sin(),1 - 2 * s **2,2 * s * (1 - s ** 2) ** rational(1,2),Omega.cos(),Omega.sin()])
		# Position vector in space.
		r = dot(Rxq,ro)
		return -r[self.__axis]
	def solve_last(self,init,t):
		from copy import deepcopy
		retval = deepcopy(init)
		retval['dlam'] = retval['dlam'] + 1. / (retval['Lam'] ** 3) * t
		return retval
	def numerical(self,tf,n):
		"""
		Return a list of n state vectors starting from t = 0 to t = tf assuming as initial conditions
		the state s0 (list of 6 numbers, rectangular position + velocity).
		"""
		from pyranha.Celmec import mdelaunay2oe, oe2s
		from scipy.integrate import odeint
		from numpy import linspace
		init_list = self.init_list[0]
		s0 = oe2s(mdelaunay2oe([init_list['Lam'],init_list['P'],init_list['Q'],init_list['dlam'] + init_list['lam0'],init_list['p'],init_list['q']]))
		def dyn(s,_):
			from numpy.linalg import norm
			axis = self.__axis
			r = norm(s[:3])
			r3 = r * r * r
			retval = [s[3],s[4],s[5],-s[0] / r3,-s[1] / r3,-s[2] / r3]
			retval[3 + axis] = retval[3 + axis] + self.__thrust
			return retval
		return odeint(dyn,s0,linspace(0,tf,n))
