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

from detail import __check_uniform_type

def r_a(e,M):
	"""
	Calculate the elliptic expansion of r/a in terms of e and M.
	"""
	__check_uniform_type(e,M)
	try:
		return e.r_a(e,M)
	except AttributeError:
		raise AttributeError, "The series type '" + str(type(e)) +  "' does not offer an r_a method."

def a_r(e,M):
	"""
	Calculate the elliptic expansion of a/r in terms of e and M.
	"""
	__check_uniform_type(e,M)
	try:
		return e.a_r(e,M)
	except AttributeError:
		raise AttributeError, "The series type '" + str(type(e)) +  "' does not offer an a_r method."

def cos_f(e,M):
	"""
	Calculate the elliptic expansion of cos(f) in terms of e and M.
	"""
	__check_uniform_type(e,M)
	try:
		return e.cos_f(e,M)
	except AttributeError:
		raise AttributeError, "The series type '" + str(type(e)) +  "' does not offer a cos_f method."

def sin_f(e,M):
	"""
	Calculate the elliptic expansion of sin(f) in terms of e and M.
	"""
	__check_uniform_type(e,M)
	try:
		return e.sin_f(e,M)
	except AttributeError:
		raise AttributeError, "The series type '" + str(type(e)) +  "' does not offer a sin_f method."

def cos_E(e,M):
	"""
	Calculate the elliptic expansion of cos(E) in terms of e and M.
	"""
	__check_uniform_type(e,M)
	try:
		return e.cos_E(e,M)
	except AttributeError:
		raise AttributeError, "The series type '" + str(type(e)) +  "' does not offer a cos_E method."

def sin_E(e,M):
	"""
	Calculate the elliptic expansion of sin(E) in terms of e and M.
	"""
	__check_uniform_type(e,M)
	try:
		return e.sin_E(e,M)
	except AttributeError:
		raise AttributeError, "The series type '" + str(type(e)) +  "' does not offer a sin_E method."

def EE(e,M):
	"""
	Calculate the elliptic expansion of E in terms of e and M.
	"""
	__check_uniform_type(e,M)
	try:
		return e.EE(e,M)
	except AttributeError:
		raise AttributeError, "The series type '" + str(type(e)) +  "' does not offer an E method."

def eipE(e,M,p = 1):
	"""
	Calculate the elliptic expansion of exp[i*p*E] in terms of e and M.
	"""
	__check_uniform_type(e,M)
	try:
		return e.eipE(e,M,p)
	except AttributeError:
		raise AttributeError, "The series type '" + str(type(e)) +  "' does not offer an eipE method."

def cos_psi():
	"""
	Calculate the expansion of cos(psi) in terms of the classical orbital elements:
	e, ep, l, lp, s, sp, v, vp, O, Op,
	where l is the mean longitude, s = sin(i/2), v the longitude of pericentre, O the longitude of ascending node.
	psi will be the angular separation of the two secondary bodies, orbiting a massive primary body, whom the orbital
	elements listed above refer to.
	"""
	from pyranha.Core import psyms, psym
	from pyranha import ds
	from pyranha.Math import cos, sin, root
	for i in "e M O s l v ep Mp Op sp lp vp".split():
		try:
			exec("%s = ds(psyms[\"%s\"])" % (i,i))
		except UserWarning:
			exec("%s = ds(psym(\"%s\"))" % (i,i))
	cos_of = cos(v-O)*cos_f(e,M).sub(psyms["M"],l-v)-sin(v-O)*sin_f(e,M).sub(psyms["M"],l-v)
	sin_of = sin(v-O)*cos_f(e,M).sub(psyms["M"],l-v)+cos(v-O)*sin_f(e,M).sub(psyms["M"],l-v)
	x_r = cos(O)*cos_of-sin(O)*sin_of*(1-2*s**2)
	y_r = sin(O)*cos_of+cos(O)*sin_of*(1-2*s**2)
	z_r = sin_of*2*s*root(2,1-s**2)
	cos_opfp = cos(vp-Op)*cos_f(ep,Mp).sub(psyms["Mp"],lp-vp)-sin(vp-Op)*sin_f(ep,Mp).sub(psyms["Mp"],lp-vp)
	sin_opfp = sin(vp-Op)*cos_f(ep,Mp).sub(psyms["Mp"],lp-vp)+cos(vp-Op)*sin_f(ep,Mp).sub(psyms["Mp"],lp-vp)
	x_rp = cos(Op)*cos_opfp-sin(Op)*sin_opfp*(1-2*sp**2)
	y_rp = sin(Op)*cos_opfp+cos(Op)*sin_opfp*(1-2*sp**2)
	z_rp = sin_opfp*2*sp*root(2,1-sp**2)
	cp = x_r*x_rp+y_r*y_rp+z_r*z_rp
	return cp

def vsop_to_dps(input_filename):
	import re
	arguments = "[poly_arg]\nname=t\ntime_eval=0;1\n"
	arguments += "[trig_arg]\nname=l_me\ntime_eval=4.40260884240;26087.9031415742\n"
	arguments += "[trig_arg]\nname=l_v\ntime_eval=3.17614669689;10213.2855462110\n"
	arguments += "[trig_arg]\nname=l_e\ntime_eval=1.75347045953;6283.0758499914\n"
	arguments += "[trig_arg]\nname=l_ma\ntime_eval=6.20347611291;3340.6124266998\n"
	arguments += "[trig_arg]\nname=l_j\ntime_eval=0.59954649739;529.6909650946\n"
	arguments += "[trig_arg]\nname=l_s\ntime_eval=0.87401675650;213.2990954380\n"
	arguments += "[trig_arg]\nname=l_u\ntime_eval=5.48129387159;74.7815985673\n"
	arguments += "[trig_arg]\nname=l_n\ntime_eval=5.31188628676;38.1330356378\n"
	arguments += "[trig_arg]\nname=D\ntime_eval=5.19846674103;77713.7714681205\n"
	arguments += "[trig_arg]\nname=F\ntime_eval=1.62790523337;84334.6615813083\n"
	arguments += "[trig_arg]\nname=l\ntime_eval=2.35555589827;83286.9142695536\n"
	arguments += "[trig_arg]\nname=Lm\ntime_eval=3.81034454697;83997.0911355954\n"
	arguments += "[terms]\n"
	l = open(input_filename).read().strip().splitlines()
	retval = []
	alpha = None
	t_rng = tuple(range(0,36,3))
	for i in l:
		tmp_header = re.match(".*\*T\*\*(\d).*",i)
		# If we found an header, read alpha and if necessary append new retval.
		if tmp_header:
			alpha = tmp_header.group(1)
			if int(alpha) == 0:
				retval.append("# "+i+"\n")
				retval[-1] += arguments
		else:
			# Do something only if alpha is not None.
			if alpha != None:
				# The multipliers are contained into the columns ranging from 10 to 10 + 12*3.
				trig_array = ";".join([i[10:46][n:n+3].strip() for n in t_rng])
				# S and K are the first two double values in the record.
				S, K = i[46:].strip().split()[0:2]
				retval[-1] += S + "!" + alpha + "|" + trig_array + ";s\n"
				retval[-1] += K + "!" + alpha + "|" + trig_array + ";c\n"
	return retval

def delaunay2oe(d_elements = None, degree = 10, t = None):
	"""
	Returns a list containing the expressions of classical orbital elements in terms of series of modified Delaunay variables.
	If d_elements is not None, all other arguments are ignored and d_elements is assumed to be a list of series representing
	the following modified Delaunay variables:

	- Lambda,
	- P,
	- Q,
	- lambda,
	- p,
	- q.

	The truncation of the series expansions needed for the transformation from Delaunay variables to classical orbital elements is
	established by the active truncator.

	If d_elements is None, then the series representing the Delaunay variables are created internally. The series will be of type t
	or qqps if t is None. Series truncation will be degree-based, with the maximum degree for the Delaunay variable P being the input
	parameter 'degree'. Note that for small eccentricities and inclinations, P is proportional to e ** 2.

	In all cases, the return value of this function will be a list representing the series expansions of the following classical
	orbital elements in terms of modifed Delaunay variables:

	- a, semi-major axis,
	- e, eccentricity,
	- sin(i/2), sine of half inclination (here referred to as 's'),
	- omega, argument of pericentre,
	- Omega, longitude of the ascending node,
	- M, mean anomaly,

	Additionally, the resulting series will contain also as symbolic variable the integral constant 2 as "two", which
	is kept in literal form since its square root appears in the expressions involving the Delaunay elements.

	Finally, note that the gravitational parameter in this transformation (e.g.,G * m, G * (m0 + m1), etc.) is
	set equal to unity.
	"""
	from pyranha.Math import root
	from pyranha.Core import psym
	if d_elements is None:
		from pyranha.Qqps import qqps
		import pyranha.Truncators
		if degree <= 0:
			raise ValueError('Truncation degree must be a positive value.')
		# st is the series type.
		if t is None:
			st = qqps
		else:
			st = t
		# Reset all the truncators.
		pyranha.Truncators.unset()
		# Create the series representing the Delaunay elements.
		Lambda = st(psym('Lambda'))
		P = st(psym('P'))
		Q = st(psym('Q'))
		lambda_ = st(psym('lambda'))
		p = st(psym('p'))
		q = st(psym('q'))
		# Now set the truncator for the series expansions.
		pyranha.Truncators.degree.set('P',degree)
	else:
		try:
			iter(d_elements)
		except TypeError:
			raise TypeError('Please provide a list of series as input parameter.')
		if len(d_elements) != 6:
			raise ValueError('The list of Delaunay variables must contain 6 elements.')
		if [type(d_elements[0])] * 6 != [type(e) for e in d_elements]:
			raise TypeError('The series representing Delaunay variables must be all of the same type.')
		st = type(d_elements[0])
		Lambda = d_elements[0]
		P = d_elements[1]
		Q = d_elements[2]
		lambda_ = d_elements[3]
		p = d_elements[4]
		q = d_elements[5]
	retval = []
	# Let's declare the "two" symbol and construct a series from it.
	two = st(psym('two',2))
	# Let's start with the semi-major axis, easy.
	a = Lambda ** 2
	retval.append(a)
	# Calculate series expansion for e.
	e = root(2, two * P * Lambda ** -1 - P ** 2 * Lambda ** -2)
	retval.append(e)
	# Now let's deal with sin(i/2) (aka 's').
	s = root(2,Q * (two * Lambda) ** -1) * root(-2, 1 - P * Lambda ** -1)
	retval.append(s)
	# Finally, let's deal with the angle variables.
	retval.append(q - p) # omega
	retval.append(-q) # Omega
	retval.append(lambda_ + p) # M
	return retval
