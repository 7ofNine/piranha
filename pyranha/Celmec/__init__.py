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

def mdelaunay2oe(md_elements):
	"""
	Transform modified Delaunay elements into classical orbital elements.

	md_elements is assumed to be a list of the following modified Delaunay variables:

	- Lambda,
	- P,
	- Q,
	- lambda,
	- p,
	- q.

	In case md_elements is a list of series, the truncation of the series expansions needed for the transformation from
	modified Delaunay variables to classical orbital elements is established by the active truncator.

	The return value of this function will be a list containing the following classical orbital element equivalent to the
	input modifed Delaunay variables:

	- a, semi-major axis,
	- e, eccentricity,
	- sin(i/2), sine of half inclination (often known as 's'),
	- omega, argument of pericentre,
	- Omega, longitude of the ascending node,
	- M, mean anomaly,

	When working with series, the resulting series will contain also as symbolic variable the integral constant 2 as "two", which
	is kept in literal form since its square root appears in the expressions involving the Delaunay elements.

	Finally, note that the gravitational parameter in this transformation (e.g.,G * m, or G * (m0 + m1), etc.) is
	set equal to unity.
	"""
	from pyranha.Math import root
	from pyranha.Core import psym, is_iteratable
	from pyranha import manipulators
	if not is_iteratable(md_elements):
		raise TypeError('Please provide a list of modified Delaunay arguments as input parameter.')
	if len(md_elements) != 6:
		raise ValueError('The list of modified Delaunay variables must contain 6 elements.')
	if [type(md_elements[0])] * 6 != [type(e) for e in md_elements]:
		raise TypeError('The modified Delaunay variables must be all of the same type.')
	t = type(md_elements[0])
	Lam, P, Q, lam, p, q = md_elements
	retval = []
	if t in manipulators:
		# Let's declare the "two" symbol and construct a series from it.
		two = t(psym('two',2))
	else:
		two = 2.
	# Let's start with the semi-major axis, easy.
	a = Lam ** 2
	retval.append(a)
	# Calculate series expansion for e.
	e = root(2, two * P * Lam ** -1 - P ** 2 * Lam ** -2)
	retval.append(e)
	# Now let's deal with sin(i/2) (aka 's').
	s = root(2,Q * (two * Lam) ** -1) * root(-2, 1 - P * Lam ** -1)
	retval.append(s)
	# Finally, let's deal with the angle variables.
	retval.append(q - p) # omega
	retval.append(-q) # Omega
	retval.append(lam + p) # M
	return retval

def oe2mdelaunay(oe):
	"""
	Transform classical orbital elements into modified Delaunay elements.

	oe is assumed to be a list of the following classical orbital elements:

	- a, semi-major axis,
	- e, eccentricity,
	- sin(i/2), sine of half inclination (often known as 's'),
	- omega, argument of pericentre,
	- Omega, longitude of the ascending node,
	- M, mean anomaly.

	In case oe is a list of series, the truncation of the series expansions needed for the transformation from
	classical orbital elements to modified Delaunay variables is established by the active truncator.

	The return value of this function will be a list containing the following modifed Delaunay variables equivalent to the
	input classical orbital element:

	- Lambda,
	- P,
	- Q,
	- lambda,
	- p,
	- q.

	Finally, note that the gravitational parameter in this transformation (e.g.,G * m, or G * (m0 + m1), etc.) is
	set equal to unity.
	"""
	from pyranha.Math import root
	from pyranha.Core import is_iteratable
	if not is_iteratable(oe):
		raise TypeError('Please provide a list of classical orbital elements as input parameter.')
	if len(oe) != 6:
		raise ValueError('The list of classical orbital elements must contain 6 elements.')
	if [type(oe[0])] * 6 != [type(e) for e in oe]:
		raise TypeError('The classical orbital elements must be all of the same type.')
	t = type(oe[0])
	a, e, s, omega, Omega, M = oe
	retval = []
	# Let's start with Lambda.
	Lam = root(2,a)
	retval.append(Lam)
	# P.
	P = Lam * (1 - root(2,1 - e * e))
	retval.append(P)
	# Q.
	Q = 2 * Lam * root(2,1 - e * e) * (s * s)
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
	from pyranha import manipulators
	if not type(s1) in manipulators or not type(s2) in manipulators:
		raise TypeError('s1 and s2 must be series.')
	if len(p_list) != len(q_list):
		raise ValueError('The list of names of momenta and coordinates must contain the same number of elements.')
	if not all(isinstance(n,str) for n in p_list) or not all(isinstance(n,str) for n in q_list):
		raise TypeError('The list of names of momenta and coordinates must contain only string elements.')
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

def lieS(eps,chi,arg,p_list,q_list):
	"""
	Lie series of generator chi developed in powers of the small quantity eps, using p_list and q_list
	as lists of names of the canonical momenta and coordinates. The limit of the series expansion is taken
	from the active truncator.
	"""
	from copy import copy
	n = eps.psi()
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

def lieH(H,eps,mom,coord,he_solver):
	"""
	Generator of chain of canonical transformations based on Lie series. Input parameters:

	- H:         series representing the full Hamiltonian
	- eps:       name of the symbol representing the small quantity
	- mom:       list of names of the momenta symbols
	- coord:     list of names of the coordinates symbols
	- limit:     desired order of quantity eps for the remainders of Lie series expansions
	- he_solver: solver for the homological equation

	he_solver must be a function (i.e., a callable) that accepts two input parameters, R and H0. R is the
	remainder of the Hamiltonian at the current order and H0 is the integrable hamiltonian. The return value
	will be the generating function of the Lie series transformation used for the next step of the iteration.
	"""
	from copy import copy
	from math import ceil
	from pyranha.Core import is_iteratable, psym
	from pyranha.Celmec import lieS
	from pyranha import manipulators
	# Check input parameters.
	if not type(H) in manipulators:
		raise TypeError('H must be a series.')
	if not isinstance(eps,str):
		raise TypeError('eps must be a string.')
	if not is_iteratable(mom) or not is_iteratable(coord):
		raise TypeError('mom and coord must be lists of strings.')
	if not all(isinstance(n,str) for n in mom) or not all(isinstance(n,str) for n in coord):
		raise TypeError('mom and coord must be lists of strings.')
	__eps = copy(eps)
	__mom = copy(mom)
	__coord = copy(coord)
	__H = copy(H)
	t = type(__H)
	s_eps = t(psym(__eps))
	n_it = s_eps.psi()
	# Integrable Hamiltonian.
	H0 = H.filtered([lambda t:(t[0] * t[1]).order(__eps) == 0] * len(H.arguments))
	# Iteration starts.
	for i in range(0,n_it - 1):
		# Residual of order i + 1.
		R = __H.filtered([lambda term:(term[0] * term[1]).order(__eps) == i + 1] * len(__H.arguments)).sub(__eps,t(1))
		chi = he_solver(R,H0)
		__H = lieS(s_eps ** (i + 1),chi,__H,__mom,__coord)
		yield copy(__H),copy(chi)

def orbitalR(angles):
	"""
	Return the rotation matrix from the orbital plane to the three-dimensional reference plane in which
	a keplerian orbit is emebedded. The rotation angles are the classical orbital elements omega, i, Omega,
	and must be passed as a list.

	If the input list contains 3 elements, they will be interpreted as the angles omega, i, Omega. If the input list
	contains 6 elements, they will be interpreted as [cos(om),sin(om),cos(i),sin(i),cos(Om),sin(Om)].
	"""
	from numpy import array
	from pyranha.Math import cos, sin
	from pyranha.Core import is_iteratable
	if not is_iteratable(angles):
		raise TypeError('Input parameter is not a list.')
	l = len(angles)
	if l != 3 and l != 6:
		raise ValueError('Input list must contain either 3 or 6 elements.')
	if l == 3:
		omega, i, Omega = angles
		cos_omega, sin_omega = cos(omega), sin(omega)
		cos_i, sin_i = cos(i), sin(i)
		cos_Omega, sin_Omega = cos(Omega), sin(Omega)
	else:
		cos_omega, sin_omega, cos_i, sin_i, cos_Omega, sin_Omega = angles
	return array([ \
		[cos_Omega * cos_omega - sin_Omega * cos_i * sin_omega, \
		-cos_Omega * sin_omega - sin_Omega * cos_i * cos_omega, \
		sin_Omega * sin_i], \
		[sin_Omega * cos_omega + cos_Omega * cos_i * sin_omega, \
		-sin_Omega * sin_omega + cos_Omega * cos_i * cos_omega, \
		-cos_Omega * sin_i], \
		[sin_i * sin_omega, \
		sin_i * cos_omega, \
		cos_i] \
	])
