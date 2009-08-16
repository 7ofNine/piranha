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
	or ds if t is None. Series truncation will be degree-based, with the maximum degree for the Delaunay variable P being the input
	parameter 'degree'. Note that for small eccentricities and inclinations, P is proportional to e ** 2.

	In all cases, the return value of this function will be a list representing the series expansions of the following classical
	orbital elements in terms of modifed Delaunay variables:

	- a, semi-major axis,
	- e, eccentricity,
	- sin(i/2), sine of half inclination (often known as 's'),
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
		from pyranha import ds
		from pyranha.Truncators import truncators
		if degree <= 0:
			raise ValueError('Truncation degree must be a positive value.')
		# st is the series type.
		if t is None:
			st = ds
		else:
			st = t
		# Reset all the truncators.
		truncators.unset()
		# Create the series representing the Delaunay elements.
		Lambda = st(psym('Lambda'))
		P = st(psym('P'))
		Q = st(psym('Q'))
		lambda_ = st(psym('lambda'))
		p = st(psym('p'))
		q = st(psym('q'))
		# Now set the truncator for the series expansions.
		truncators.degree.set('P',degree)
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

def oe2delaunay(oe = None, degree = 10, t = None):
	"""
	Returns a list containing the expressions of modified Delaunay variables in terms of classical orbital elements.
	If oe is not None, all other arguments are ignored and oe is assumed to be a list of series representing
	the following classical orbital elements:

	- a, semi-major axis,
	- e, eccentricity,
	- sin(i/2), sine of half inclination (often known as 's'),
	- omega, argument of pericentre,
	- Omega, longitude of the ascending node,
	- M, mean anomaly,

	The truncation of the series expansions needed for the transformation from classical orbital elements to modified Delaunay variables
	is established by the active truncator.

	If oe is None, then the series representing the classical orbital elements are created internally. The series will be of type t
	or ds if t is None. Series truncation will be degree-based, with the maximum degree for the eccentricity being the input
	parameter 'degree'.

	In all cases, the return value of this function will be a list representing the series expansions of the following modified Delaunay
	variables in terms of classical orbital elements:

	- Lambda,
	- P,
	- Q,
	- lambda,
	- p,
	- q.

	Note that the gravitational parameter in this transformation (e.g.,G * m, G * (m0 + m1), etc.) is set equal to unity.
	"""
	from pyranha.Math import root
	if oe is None:
		from pyranha.Core import psym
		from pyranha import ds
		from pyranha.Truncators import truncators
		if degree <= 0:
			raise ValueError('Truncation degree must be a positive value.')
		# st is the series type.
		if t is None:
			st = ds
		else:
			st = t
		# Reset all the truncators.
		truncators.unset()
		# Create the series representing the classical orbital elements.
		a = st(psym('a'))
		e = st(psym('e'))
		s = st(psym('s'))
		omega = st(psym('omega'))
		Omega = st(psym('Omega'))
		M = st(psym('M'))
		# Now set the truncator for the series expansions.
		truncators.degree.set('e',degree)
	else:
		try:
			iter(oe)
		except TypeError:
			raise TypeError('Please provide a list of series as input parameter.')
		if len(oe) != 6:
			raise ValueError('The list of classical orbital elements must contain 6 elements.')
		if [type(oe[0])] * 6 != [type(e) for e in oe]:
			raise TypeError('The series representing classical orbital elements must be all of the same type.')
		st = type(oe[0])
		a = oe[0]
		e = oe[1]
		s = oe[2]
		omega = oe[3]
		Omega = oe[4]
		M = oe[5]
	retval = []
	# Let's start with Lambda.
	Lambda = root(2,a)
	retval.append(Lambda)
	# P.
	P = Lambda * (1 - root(2,1 - e * e))
	retval.append(P)
	# Q.
	Q = 2 * Lambda * root(2,1 - e * e) * (s * s)
	retval.append(Q)
	# Finally, let's deal with the angle variables.
	retval.append(M + omega + Omega) # lambda
	retval.append(-omega -Omega) # p
	retval.append(-Omega) # q
	return retval

def poisson_bra(s1,s2,p_list,q_list):
	"""
	Calculate the Poisson bracket {s1,s2}, with series s1 and s2 are interpreted as
	functions of generalised momenta whose names are listed in p_list and generalised
	coordinates whose names are listed in q_list.
	"""
	from pyranha.Math import partial
	if len(p_list) != len(q_list):
		raise ValueError('The list of names of momenta and coordinates must contain the same number of elements.')
	for i in p_list:
		if not isinstance(i,str):
			raise TypeError('The list of names of momenta variables must contain only string elements.')
	for i in q_list:
		if not isinstance(i,str):
			raise TypeError('The list of names of coordinate variables must contain only string elements.')
	l = [partial(s1,q_list[i]) * partial(s2,p_list[i]) - partial(s1,p_list[i]) * partial(s2,q_list[i]) for i in range(0,len(p_list))]
	return sum(l)

def is_canonical(new_p,new_q,p_list,q_list):
	"""
	Test whether the transformation new_p = new_p(p_list) and new_q = new_q(q_list)
	is canonical using the Poisson brackets test.

	new_p and new_q must be lists of series, while p_list and q_list must be lists of
	variable names.
	"""
	n = len(new_p)
	if n != len(new_q) or n != len(p_list) or n != len(q_list):
		raise ValueError('The lengths of the input lists must be consistent.')
	for i in range(0,n):
		for j in range(0,n):
			if poisson_bra(new_p[i],new_p[j],p_list,q_list) != 0:
				return False
			if poisson_bra(new_q[i],new_q[j],p_list,q_list) != 0:
				return False
			if poisson_bra(new_q[i],new_p[j],p_list,q_list) != int(i == j):
				return False
	return True

def lieL(gen,arg,p_list,q_list,n = 1):
	"""
	Lie derivative of order n on argument arg with generator gen, using p_list and q_list
	as lists of names of the canonical momenta and coordinates.
	"""
	if n < 0:
		raise ValueError('The order of the Lie operator must be non-negative.')
	if n == 0:
		return arg
	retval = poisson_bra(arg,gen,p_list,q_list)
	for _ in range(1,n):
		retval = poisson_bra(retval,gen,p_list,q_list)
	return retval

def lieS(eps,chi,arg,p_list,q_list,order = None):
	"""
	Lie series of generator chi developed in powers of the small quantity eps, using p_list and q_list
	as lists of names of the canonical momenta and coordinates.

	If order is None, then the limit of the series is taken from the active truncator. Otherwise, the series limit
	is given by order itself (which must be a non-negative integer) and which represents the power of eps from which the
	remainder of the Lie series starts.
	"""
	from copy import copy
	if order is None:
		n = eps.psi()
	else:
		n = order
	if n == 0:
		return type(eps)()
	tmp = copy(arg)
	retval = copy(arg)
	for i in range(1,n):
		tmp = lieL(chi,tmp,p_list,q_list)
		tmp *= eps
		tmp /= i
		retval += tmp
	return retval

def orbitalR(angles, t = None):
	"""
	Return the rotation matrix from the orbital plane to the three-dimensional reference plane in which
	a keplerian orbit is emebedded. The rotation angles are the classical orbital elements omega, i, Omega,
	and must be passed as a list.
	"""
	from numpy import array
	from pyranha.Math import cos, sin
	omega, i, Omega = angles
	return array([ \
		[cos(Omega) * cos(omega) - sin(Omega) * cos(i) * sin(omega), \
		-cos(Omega) * sin(omega) - sin(Omega) * cos(i) * cos(omega), \
		sin(Omega) * sin(i)], \
		[sin(Omega) * cos(omega) + cos(Omega) * cos(i) * sin(omega), \
		-sin(Omega) * sin(omega) + cos(Omega) * cos(i) * cos(omega), \
		-cos(Omega) * sin(i)], \
		[sin(i) * sin(omega), \
		sin(i) * cos(omega), \
		cos(i)] \
	])
