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

import pyranha as __pyranha

def r_a(e,M):
	"""
	Calculate the elliptic expansion of r/a in terms of e and M.
	"""
	try:
		return __pyranha.ds.r_a(e,M)
	except AttributeError:
		raise AttributeError, "The series type '" + str(__pyranha.ds) +  "' does not offer an r_a method."

def cos_f(e,M):
	"""
	Calculate the elliptic expansion of cos(f) in terms of e and M.
	"""
	try:
		return __pyranha.ds.cos_f(e,M)
	except AttributeError:
		raise AttributeError, "The series type '" + str(__pyranha.ds) +  "' does not offer a cos_f method."

def sin_f(e,M):
	"""
	Calculate the elliptic expansion of sin(f) in terms of e and M.
	"""
	try:
		return __pyranha.ds.sin_f(e,M)
	except AttributeError:
		raise AttributeError, "The series type '" + str(__pyranha.ds) +  "' does not offer a sin_f method."

def cos_E(e,M):
	"""
	Calculate the elliptic expansion of cos(E) in terms of e and M.
	"""
	try:
		return __pyranha.ds.cos_E(e,M)
	except AttributeError:
		raise AttributeError, "The series type '" + str(__pyranha.ds) +  "' does not offer a cos_E method."

def sin_E(e,M):
	"""
	Calculate the elliptic expansion of sin(E) in terms of e and M.
	"""
	try:
		return __pyranha.ds.sin_E(e,M)
	except AttributeError:
		raise AttributeError, "The series type '" + str(__pyranha.ds) +  "' does not offer a sin_E method."

def EE(e,M):
	"""
	Calculate the elliptic expansion of E in terms of e and M.
	"""
	try:
		return __pyranha.ds.E(e,M)
	except AttributeError:
		raise AttributeError, "The series type '" + str(__pyranha.ds) +  "' does not offer an E method."

def cos_psi():
	"""
	Calculate the expansion of cos(psi) in terms of the classical orbital elements:
	e, e', l, l', s, s', v, v', O, O',
	where l is the mean longitude, s = sin(i/2), v the longitude of pericentre, O the longitude of ascending node.
	psi will be the angular separation of the two secondary bodies, orbiting a massive primary body, whom the orbital
	elements listed above refer to.
	"""
	from pyranha.Core import psym as psym
	from pyranha import ds as ds
	from pyranha.Math import cos as cos, sin as sin, root as root
	e,M,O,s,l,v = map(psym,'eMOslv')
	ep,Mp,Op,sp,lp,vp = map(psym,["e'","M'","O'","s'","l'","v'"])
	cos_of = cos(ds(v)-ds(O))*cos_f(e,M).sub(M,ds(l)-ds(v))-sin(ds(v)-ds(O))*sin_f(e,M).sub(M,ds(l)-ds(v))
	sin_of = sin(ds(v)-ds(O))*cos_f(e,M).sub(M,ds(l)-ds(v))+cos(ds(v)-ds(O))*sin_f(e,M).sub(M,ds(l)-ds(v))
	x_r = cos(ds(O))*cos_of-sin(ds(O))*sin_of*(1-2*ds(s)**2)
	y_r = sin(ds(O))*cos_of+cos(ds(O))*sin_of*(1-2*ds(s)**2)
	z_r = sin_of*2*ds(s)*root(2,1-ds(s)**2)
	cos_opfp = cos(ds(vp)-ds(Op))*cos_f(ep,Mp).sub(Mp,ds(lp)-ds(vp))-sin(ds(vp)-ds(Op))*sin_f(ep,Mp).sub(Mp,ds(lp)-ds(vp))
	sin_opfp = sin(ds(vp)-ds(Op))*cos_f(ep,Mp).sub(Mp,ds(lp)-ds(vp))+cos(ds(vp)-ds(Op))*sin_f(ep,Mp).sub(Mp,ds(lp)-ds(vp))
	x_rp = cos(ds(Op))*cos_opfp-sin(ds(Op))*sin_opfp*(1-2*ds(sp)**2)
	y_rp = sin(ds(Op))*cos_opfp+cos(ds(Op))*sin_opfp*(1-2*ds(sp)**2)
	z_rp = sin_opfp*2*ds(sp)*root(2,1-ds(sp)**2)
	cp = x_r*x_rp+y_r*y_rp+z_r*z_rp
	return cp
