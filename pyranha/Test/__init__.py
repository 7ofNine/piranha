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
				# Here we should wrt standard IEEE floating point format, the factorial is small enough
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
