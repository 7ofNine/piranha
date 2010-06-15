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

def mdelaunay2poincare(md):
	"""
	Transform modified Delaunay orbital elements into Poincare' variables.

	Return values: [Lambda,y,z,lambda,x,v]
	"""
	from copy import deepcopy
	from pyranha import manipulators
	from pyranha.Core import psym
	from pyranha.Math import root, cos, sin
	Lam, P, Q, lam, p, q = md
	t = type(Lam)
	if t in manipulators:
		# Let's declare the "two" symbol and construct a series from it.
		two = t(psym('two',2))
	else:
		two = 2.
	sqrt2P, sqrt2Q = root(2,two * P), root(2,two * Q)
	return [deepcopy(Lam),sqrt2P * cos(p),sqrt2Q * cos(q),deepcopy(lam),sqrt2P * sin(p),sqrt2Q * sin(q)]

def poincare2mdelaunay(pv):
	"""
	Transform Poincare' variables into modified Delaunay orbital elements.

	Return values: [Lambda,P,Q,lambda,p,q]
	"""
	from copy import deepcopy
	from math import sqrt, atan2
	Lam, y, z, lam, x, v = pv
	P = (x ** 2 + y ** 2) / 2.
	Q = (z ** 2 + v ** 2) / 2.
	p = atan2(x / sqrt(2. * P), y / sqrt(2. * P))
	q = atan2(v / sqrt(2. * Q), z / sqrt(2. * Q))
	return deepcopy([Lam,P,Q,lam,p,q])

def oe2s(oe):
	"""
	Transform classical orbital elements into state vector.
	"""
	# TODO: proper docs.
	from scipy.optimize import fsolve
	from math import sin, cos, sqrt, asin
	from numpy import dot, concatenate
	a, e, s, om, Om, M = oe
	n = sqrt(1. / (a ** 3))
	E = fsolve(lambda x: x - e * sin(x) - M, 0)
	i = 2 * asin(s)
	# Position and velocity in the orbital plane.
	rp = [a * (cos(E) - e), a * sqrt(1 - e ** 2) * sin(E), 0]
	vp = [-(n * a * sin(E)) / (1 - e * cos(E)), n * a * sqrt(1 - e ** 2) * cos(E) / (1 - e * cos(E)),0]
	# Rotate the orbital plane
	Rxq = orbitalR([om,i,Om])
	return concatenate((dot(Rxq,rp),dot(Rxq,vp)))

def s2oe(s):
	"""
	Transform state vector into classical orbital elements.
	"""
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
	Test whether the transformation new_p = new_p(p_list,q_list) and new_q = new_q(p_list,q_list)
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

def lieS(eps,gen,arg,p_list,q_list,limit = None):
	"""
	Lie series of generator gen developed in powers of the small quantity eps, using p_list and q_list
	as lists of names of the canonical momenta and coordinates. The limit of the series expansion,if not specified,
	is taken from the active truncator.
	"""
	from copy import deepcopy
	if limit is None:
		n = eps.psi()
	else:
		if not isinstance(limit,int):
			raise TypeError('The limit of thse Lie series must be an integer.')
		n = limit
	if n == 0:
		return type(eps)()
	tmp = deepcopy(arg)
	retval = deepcopy(arg)
	for i in range(1,n):
		tmp = lieL(gen,tmp,p_list,q_list)
		tmp *= eps
		tmp /= i
		retval += tmp
	return retval

class lie_cache(object):
	"""
	Cache of Lie derivatives. Will remember the Lie derivatives of lower orders and use them to calculate higher order derivatives.
	"""
	def __init__(self,gen,arg,p_list,q_list):
		"""
		Initialise with generator gen, argument arg, list of momenta names p_list and list of coordinates names q_list.
		"""
		from copy import deepcopy
		self.__gen = deepcopy(gen)
		self.__p_list = deepcopy(p_list)
		self.__q_list = deepcopy(q_list)
		self.__container = [deepcopy(arg)]
	def __getitem__(self,n):
		"""
		Require the Lie derivative of order n.
		"""
		if not isinstance(n,int) or n < 0:
			raise ValueError('Invalid order.')
		from copy import deepcopy
		for i in range(n - len(self.__container) + 1):
			self.__container.append(lieL(self.__gen,self.__container[-1],self.__p_list,self.__q_list))
		return deepcopy(self.__container[n])

class lie_theory(object):
	"""
	Lie theory class.
	
	This class produces a Lie theory given an Hamiltonian, a set of variable names and a list of
	solvers for the homological equations.
	
	Each class deriving from lie_theory must implement the abstract method solve_last().
	"""
	from abc import ABCMeta as __abc_meta
	from abc import abstractmethod as __abs_method
	__metaclass__ = __abc_meta
	def __init__(self,H,eps_name,p_names,q_names,he_solvers):
		"""
		Construct a Lie theory from:
		
		H -- series representing the Hamiltonian
		eps_name -- name of the variable representing the small quantity (eps)
		p_names -- list of names of variables representing the momenta
		q_names -- list of names of variables representing the coordinates
		he_solvers -- list of callable object representing the homological equation solvers
		
		The final order of the theory is given by len(he_solvers). For each order n of the theory,
		a homological equation solver is intended to take as argument the part of the current Hamiltonian
		which contains terms in eps ** n and return the generating function for the transformation.
		"""
		from copy import deepcopy
		from pyranha import manipulators
		from pyranha.Core import psym
		# Input parameters checks.
		if not type(H) in manipulators:
			raise(TypeError('Hamiltonian must be a series.'))
		self.__H = deepcopy(H)
		self.__series_type = type(H)
		if not isinstance(eps_name,str) or not all([isinstance(i,str) for i in p_names]) or not all([isinstance(i,str) for i in q_names]):
			raise(TypeError('Variable names must be strings.'))
		self.__eps_name = deepcopy(eps_name)
		self.__p_names = list(p_names)
		self.__q_names = list(q_names)
		if len(self.__p_names) != len(self.__q_names):
			raise(ValueError('The numbers of momenta and coordinates must be the same.'))
		if len(self.__p_names) == 0:
			raise(ValueError('Problem dimension must be strictly positive.'))
		self.__he_solvers = list(he_solvers)
		if not all([callable(i) for i in self.__he_solvers]):
			raise(TypeError('Homological equation solvers must all be callable objects.'))
		if len(self.__he_solvers) < 1:
			raise(ValueError('The order of the theory must be at least 1.'))
		self.__order = len(self.__he_solvers)
		self.__init_list = None
		# Run the calculations.
		self.__compute_theory()
	def __repr__(self):
		retval = ''
		retval += 'Hamiltonian:\n\t' + str(self.__H) + '\n\n'
		retval += 'Eps:\n\t' + '\'' + str(self.__eps_name) + '\'\n'
		retval += 'Momenta:\n\t' + str(self.__p_names) + '\n'
		retval += 'Coordinates:\n\t' + str(self.__q_names) + '\n'
		retval += 'Order:\n\t' + str(self.__order) + '\n'
		return retval
	def __compute_theory(self):
		from pyranha.Core import psym, integer
		self.__chi_list = []
		self.__direct = []
		self.__inverse = []
		self.__H_list = [self.__H]
		eps_series = self.__series_type(psym(self.__eps_name))
		state_series = [self.__series_type(psym(s)) for s in self.__p_names] + [self.__series_type(psym(s)) for s in self.__q_names]
		for i in range(0,self.__order):
			print('\033[1;32mCurrent perturbative order: %d of %d\033[1;m' % (i + 1,self.__order))
			print('Decomposing current Hamiltonian...')
			# TODO: fix the filtering.
			H_split = [self.__H_list[i].filtered([None, lambda t: t[1].degree(self.__eps_name) == j]).sub(self.__eps_name,self.__series_type(1)) for j in range(0,self.__order + 1)]
			print('Solving homological equation...')
			chi_n = self.__he_solvers[i](H_split[i + 1])
			print('Initialising Hamiltonian Lie caches...')
			H_lie_caches = [lie_cache(chi_n,Hn,self.__p_names,self.__q_names) for Hn in H_split]
			print('Calculating Hamiltonian...')
			H_n = self.__series_type(0)
			for j in range(0,self.__order + 1):
				tmp = self.__series_type(0)
				for k in range(0,j / (i + 1) + 1):
					tmp += H_lie_caches[j - k * (i + 1)][k] / (integer(k).factorial())
				tmp *= eps_series ** j
				H_n += tmp
			self.__H_list.append(H_n)
			self.__chi_list.append(chi_n)
			s_limit = self.__order / (i + 1) + 1
			print('Calculating ' +'\033[1;31mdirect\033[1;m ' +'transformations...')
			self.__direct.append([lieS(eps_series ** (i + 1),chi_n,x,self.__p_names,self.__q_names,limit = s_limit) for x in state_series])
			print('Calculating ' +'\033[1;31minverse\033[1;m ' +'transformations...')
			self.__inverse.append([lieS(eps_series ** (i + 1),-chi_n,x,self.__p_names,self.__q_names,limit = s_limit) for x in state_series])
			print(H_n.filtered([None, lambda t: t[1].degree(self.__eps_name) < i + 2]))
	@__abs_method
	def solve_last(self,init,t):
		"""
		Given a dictionary of initial conditions for all symbolic variables, return a dictionary with the value of
		symbolic variables at time t.
		
		init -- dictionary with the values of the symbolic variables at time t = 0
		t -- time the return dictionary refers to
		"""
		pass
	def set_init(self,init):
		"""
		Set initial conditions of the theory.
		
		init -- dictionary with the values of the symbolic variables at time t = 0
		"""
		from copy import deepcopy
		s_names = self.__p_names + self.__q_names
		# Check that init eval dictionary has all the necessary variables.
		if not all([s in init for s in s_names]):
			raise(ValueError('Initial conditions are needed for *all* coordinates and momenta.'))
		init_list = [deepcopy(init)]
		for i in range(0,self.__order):
			# Copy the original eval dictionary and transform coordinates and momenta.
			cur_init = deepcopy(init)
			for j in range(0,len(s_names)):
				cur_init[s_names[j]] = self.__inverse[i][j].eval(init_list[-1])
			init_list.append(cur_init)
		self.__init_list = init_list
	def evaluate(self,t):
		"""
		Return dictionary of symbolic variables evaluated at time t.
		
		t -- time
		"""
		from copy import deepcopy
		s_names = self.__p_names + self.__q_names
		evals = [self.solve_last(self.init_list[-1],t)]
		for d in list(reversed(self.__direct)):
			cur_eval = deepcopy(evals[-1])
			for j in range(0,len(s_names)):
				cur_eval[s_names[j]] = d[j].eval(evals[-1])
			evals.append(cur_eval)
		return evals[-1]
	@property
	def H(self):
		"""
		List of Hamiltonians for all the orders of the theory.
		"""
		return self.__H_list
	@property
	def direct(self):
		"""
		List of direct coordinate transformations for all the orders of the theory.
		"""
		return self.__direct
	@property
	def inverse(self):
		"""
		List of inverse coordinate transformations for all the orders of the theory.
		"""
		return self.__inverse
	@property
	def chi(self):
		"""
		List of generating functions for all the orders of the theory.
		"""
		return self.__chi_list
	@property
	def init_list(self):
		"""
		List of initial conditions for all the orders of the theory.
		"""
		from copy import deepcopy
		if self.__init_list is None:
			raise(ValueError('Initial conditions have not been set. Please use set_init().'))
		return deepcopy(self.__init_list)
	@property
	def series_type(self):
		"""
		Series type used in the theory.
		"""
		return self.__series_type
	@property
	def eps_name(self):
		"""
		Name of the small quantity.
		"""
		return self.__eps_name
	@property
	def p_names(self):
		"""
		List of names of the momenta.
		"""
		return self.__p_names
	@property
	def q_names(self):
		"""
		List of names of the coordinates.
		"""
		return self.__q_names
	@property
	def order(self):
		"""
		Order of the theory.
		"""
		return self.__order

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
