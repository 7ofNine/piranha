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

import unittest
import pyranha.Core

integer_numerical_types = [int, pyranha.Core.integer]

scalar_exact_numerical_types = [pyranha.Core.integer, pyranha.Core.rational]
scalar_numerical_types = integer_numerical_types + [float, pyranha.Core.rational]

complex_exact_numerical_types = []
complex_numerical_types = [complex]

exact_numerical_types = scalar_exact_numerical_types + complex_exact_numerical_types
numerical_types = scalar_numerical_types + complex_numerical_types

class double_factorial_test(unittest.TestCase):
	"""
	Exercise known relations between double factorials and factorials for integer numerical types.
	"""
	def runTest(self):
		from pyranha.Math import factorial, double_factorial
		from copy import copy
		for t in integer_numerical_types:
			self.assertEqual(double_factorial(t(0)),1)
			self.assertEqual(double_factorial(t(1)),1)
			self.assertEqual(factorial(t(0)),1)
			self.assertEqual(factorial(t(1)),1)
			for n in range(2,10):
				# Some known relations involving factorials and double factorials.
				self.assertEqual(double_factorial(t(n) - 2),double_factorial(t(n)) / n)
				self.assertEqual(factorial(t(n)),double_factorial(t(n)) * double_factorial(t(n) - 1))
				self.assertEqual(double_factorial(t(n) * 2 + 1),factorial(t(n) * 2 + 1) / double_factorial(t(n) * 2))
				self.assertEqual(double_factorial(t(n) * 2 + 1),factorial(t(n) * 2 + 1) / (t(2) ** n * factorial(t(n))))
				self.assertEqual(double_factorial(t(n) * 2 - 1),factorial(t(n) * 2 - 1) / double_factorial(t(n) * 2 - 2))
				self.assertEqual(double_factorial(t(n) * 2 - 1),factorial(t(n) * 2) / (t(2) ** n * factorial(t(n))))
				# Check that results are the same for all types.
				t_list = copy(integer_numerical_types)
				t_list.remove(t)
				for t2 in t_list:
					self.assertEqual(double_factorial(t(n)),double_factorial(t2(n)))
					self.assertEqual(factorial(t(n) - 2),factorial(t2(n) - 2))

class rf_factorial_test(unittest.TestCase):
	"""
	Exercise known relations between rising/falling factorials, factorial and choose functions.
	"""
	def runTest(self):
		from pyranha.Math import factorial, r_factorial, f_factorial, choose, cs_phase
		from numpy import linspace
		for t in numerical_types:
			for n in range(0,8):
				# Here we should be safe wrt standard IEEE floating point format, the factorial is small enough
				# to be contained within 2**52 and we are working only with integers and multiplications.
				self.assertEqual(r_factorial(t(1),n),f_factorial(t(n),n))
				if t in integer_numerical_types:
					self.assertEqual(r_factorial(t(1),n),factorial(t(n)))
				for x in linspace(t(-10),t(10),26):
					if t in exact_numerical_types:
						self.assertEqual(r_factorial(-t(x),n),f_factorial(t(x),n) * cs_phase(n))
						self.assertEqual(r_factorial(t(x),n) / factorial(n),choose(t(x) + n - 1, n))
						self.assertEqual(f_factorial(t(x),n) / factorial(n),choose(t(x), n))
						# Binomial-like theorem for rising/falling factorials.
						for y in linspace(t(-10),t(10),26):
							self.assertEqual(r_factorial(t(x) + t(y),n),sum([choose(t(n),j) * r_factorial(t(x),n - j) * r_factorial(t(y),j) for j in range(0,n + 1)]))
							self.assertEqual(f_factorial(t(x) + t(y),n),sum([choose(t(n),j) * f_factorial(t(x),n - j) * f_factorial(t(y),j) for j in range(0,n + 1)]))

class binomial_test(unittest.TestCase):
	"""
	Exercise known relations for the binomial coefficient.
	"""
	def runTest(self):
		from pyranha.Math import choose, r_factorial, cs_phase
		from numpy import linspace
		for t in numerical_types:
			for n in range(0,10):
				self.assertEqual(choose(t(n),0),choose(t(n),n))
				self.assertEqual(choose(t(n),0),1)
				self.assertEqual(choose(t(n),-1),choose(t(n),n + 1))
				self.assertEqual(sum([choose(t(n),k) for k in range(0,n + 1)]),t(2) ** n)
				if n >= 1:
					self.assertEqual(sum([choose(t(n),k) * k for k in range(1,n + 1)]),t(2) ** (n - 1) * n)
					self.assertEqual(sum([choose(t(n),i) ** 2 * i * i for i in range(0,n + 1)]),choose(t(n) * 2 - 2, n - 1) * n * n)
				self.assertEqual(sum([choose(t(n),k) ** 2 for k in range(0,n + 1)]),choose(t(n) * 2, n))
				self.assertEqual(sum([choose(t(n),i) ** 2 * i for i in range(0,n + 1)]),choose(t(n) * 2, n) / 2 * n)
				for k in range(1,n):
					self.assertEqual(sum([choose(t(n),j) * cs_phase(j) * r_factorial(t(j),k) for j in range(0,n + 1)]),0)
				for k in range(1,10):
					self.assertEqual((t(n) - 2 * k) * choose(t(n),k), t(n) * (choose(t(n) - 1, k) - choose(t(n) - 1, k - 1)))

def suite_math():
	suite = unittest.TestSuite()
	suite.addTest(double_factorial_test())
	suite.addTest(binomial_test())
	suite.addTest(rf_factorial_test())
	return suite

from pyranha import manipulators

scalar_exact_series_types = filter(lambda m: hasattr(m,'is_ring_exact') and not hasattr(m,'is_complex'),manipulators)
complex_exact_series_types = filter(lambda m: hasattr(m,'is_ring_exact') and hasattr(m,'is_complex'),manipulators)
exact_series_types = scalar_exact_series_types + complex_exact_series_types

scalar_trig_exact_series_types = filter(lambda m: hasattr(m,'is_ring_exact') and hasattr(m,'is_trig_exact') and not hasattr(m,'is_complex'),manipulators)
complex_trig_exact_series_types = filter(lambda m: hasattr(m,'is_ring_exact') and hasattr(m,'is_trig_exact') and hasattr(m,'is_complex'),manipulators)
trig_exact_series_types = scalar_trig_exact_series_types + complex_trig_exact_series_types

class series_sf_test01(unittest.TestCase):
	"""
	Exercise known relations involving special functions, part 1.
	"""
	def runTest(self):
		# TODO: integration tests.
		from pyranha.Core import psym, rational
		from pyranha.Math import cs_phase
		from pyranha.Truncators import truncators
		from detail import check_order
		truncators.unset()
		for limit in [0,1,2,3,80]:
			truncators.degree.set(limit)
			for t in exact_series_types:
				x = t(psym('x'))
				self.assertEqual(x.root(1),x)
				self.assertEqual(x ** rational(1,1),x)
				for n in range(-10,11):
					self.assert_(t().besselJ(n) == 1 or t().besselJ(n) == 0)
					self.assertEqual(x.besselJ(n),x.besselJ(-n) * cs_phase(n))
					for m in range(-10,11):
						if n >= m:
							Jnx = x.besselJ(n)
							Jnxm = x.besselJ_div_m(n,m)
							if m >= 0:
								self.assert_(check_order(x ** m * Jnxm, Jnx,limit))
							elif m < 0:
								self.assert_(check_order(Jnxm,Jnx * x ** (-m),limit))
					if n != 0:
						self.assertEqual((x ** -n).root(-n),x)
						self.assertEqual((x ** -n) ** rational(1,-n),x)

class series_sf_test02(unittest.TestCase):
	"""
	Exercise known relations involving special functions, part 2.
	"""
	def runTest(self):
		# TODO: integration tests, like orthogonality of Legendre polynomials/functions.
		from pyranha.Core import psym, rational, integer
		from pyranha.Math import cs_phase
		from pyranha.Truncators import truncators
		# Legendre polynomials.
		truncators.unset()
		for t in filter(lambda m: hasattr(m,'is_divint_exact'), exact_series_types):
			x = t(psym('x'))
			for n in range(0,50):
				print n
				self.assertEqual((-x).legendrePn(n),cs_phase(n) * x.legendrePn(n))
				self.assertEqual(t(1).legendrePn(n),1)
				# TODO: modify this when we implement substitution by complex series.
				if not hasattr(t,'real'):
					self.assertEqual(x.legendrePn(n).partial('x').sub('x',t(1)),rational(n)*(rational(n) + 1) / 2)
				if n > 0:
					self.assertEqual((n + 1) * x.legendrePn(n + 1), (2 * n + 1) * x * x.legendrePn(n) - n * x.legendrePn(n - 1))
					self.assertEqual((x ** 2 - 1) * x.legendrePn(n).partial('x'), n * x * x.legendrePn(n) - n * x.legendrePn(n - 1))
				# Satisfy Legendre's differential equation.
				self.assertEqual(((1 - x ** 2) * x.legendrePn(n).partial('x')).partial('x') + rational(n) * (n + 1) * x.legendrePn(n), 0)
				# Another differential operation.
				if n > 0:
					self.assertEqual((2 * n + 1) * x.legendrePn(n),(x.legendrePn(n + 1) - x.legendrePn(n - 1)).partial('x'))
				# Rodrigues' formula.
				self.assertEqual(x.legendrePn(n),rational(1) / (integer(2) ** n * integer(n).factorial()) * ((x ** 2 - 1) ** n).partial('x',n))
		# Associated Legendre functions.
		# TODO: Gaunt's formula and odd values of m.
		for t in filter(lambda m: hasattr(m,'is_divint_exact'), exact_series_types):
			x = t(psym('x'))
			for n in range(-20,20):
				for m in range(-n,n):
					if n >= 0 and m >= 0 and not (m & 1):
						self.assertEqual(x.legendrePnm(n,m),cs_phase(m) * (1 - x ** 2) ** (m / 2) * x.legendrePn(n).partial('x',m))
						self.assertEqual(x.legendrePnm(n,m),rational(cs_phase(m)) / (integer(2) ** n * integer(n).factorial()) * (1 - x ** 2) ** (m / 2) * ((x ** 2 - 1) ** n).partial('x',n + m))
					if not (m & 1):
						self.assertEqual(x.legendrePnm(n,-m),cs_phase(m) * rational(integer(n - m).factorial(),integer(n + m).factorial()) * x.legendrePnm(n,m))

class series_trig_test(unittest.TestCase):
	"""
	Exercise known relations involving trigonometric functions.
	"""
	def runTest(self):
		# TODO: legendre tests with sine/cosine as arguments.
		from pyranha.Core import psym, integer, rational
		from pyranha.Math import choose, einpi2, cs_phase
		from pyranha.Truncators import truncators
		truncators.unset()
		for limit in [1,2,3,80]:
			truncators.degree.set(limit)
			for t in scalar_trig_exact_series_types:
				x = t(psym('x'))
				self.assertEqual(x.sin() * x.sin() + x.cos() * x.cos(), 1)
				# Double angle formulas.
				self.assertEqual((2 * x).sin(), 2 * x.sin() * x.cos())
				self.assertEqual((2 * x).cos(), x.cos() ** 2 - x.sin() ** 2)
				self.assertEqual((2 * x).cos(), 2 * x.cos() ** 2 - 1)
				self.assertEqual((2 * x).cos(), 1 - 2 * x.sin() ** 2)
				# Triple angle formulas.
				self.assertEqual((3 * x).sin(), 3 * x.sin() - 4 * x.sin() ** 3)
				self.assertEqual((3 * x).cos(), 4 * x.cos() ** 3 - 3 * x.cos())
				for n in range(0,21):
					# Sine/cosine of multiple angles.
					self.assertEqual((n * x).sin(), sum([choose(integer(n),k) * x.cos() ** k * x.sin() ** (n - k) * einpi2(n - k).imag for k in range(0,n + 1)]))
					self.assertEqual((n * x).cos(), sum([choose(integer(n),k) * x.cos() ** k * x.sin() ** (n - k) * einpi2(n - k).real for k in range(0,n + 1)]))
					# Power-reduction formulas.
					if n % 2:
						self.assertEqual(x.cos() ** n, rational(2) / (rational(2) ** n) * sum([choose(rational(n),k) * ((n - 2 * k) * x).cos() for k in range(0,(n - 1) / 2 + 1)]))
						self.assertEqual(x.sin() ** n, rational(2) / (rational(2) ** n) * sum([cs_phase((n - 1) / 2 - k) * choose(rational(n),k) * ((n - 2 * k) * x).sin() for k in range(0,(n - 1) / 2 + 1)]))
					else:
						self.assertEqual(x.cos() ** n, rational(1) / (rational(2) ** n) * choose(rational(n), (n / 2)) + rational(2) / (rational(2) ** n) * sum([choose(rational(n),k) * ((n - 2 * k) * x).cos() for k in range(0,n / 2)]))
						self.assertEqual(x.sin() ** n, rational(1) / (rational(2) ** n) * choose(rational(n), (n / 2)) + rational(2) / (rational(2) ** n) * sum([cs_phase(n / 2 - k) * choose(rational(n),k) * ((n - 2 * k) * x).cos() for k in range(0,n / 2)]))

def suite_series():
	suite = unittest.TestSuite()
	suite.addTest(series_sf_test01())
	suite.addTest(series_sf_test02())
	suite.addTest(series_trig_test())
	return suite

def run_full_suite():
	alltests = unittest.TestSuite([suite_series(), suite_math()])
	unittest.TextTestRunner(verbosity = 2).run(alltests)
