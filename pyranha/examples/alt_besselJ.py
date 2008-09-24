def alt_besselJ(n,x):
	"""
	Python implementation of the power series expansion for Bessel functions of the first kind
	of order n and argument x.
	"""
	# Import the Pyranha module and the "factorial" function from its Math submodule.
	import pyranha
	from pyranha.Math import factorial
	# First we check that n is an integer.
	if type(n) != int:
		raise TypeError("n must be integer.")
	# Next we check that n has a sane value.
	if n < 0:
		raise ValueError("n must be non-negative.")
	# Last, check that the type of x belongs to the available Pyranha manipulators.
	if type(x) not in pyranha.manipulators_type_tuple:
		raise TypeError("x is not a Pyranha series type.")
	# Now let's proceed to the actual expansion. Formula is available on Wikipedia, for reference:
	# http://en.wikipedia.org/wiki/Bessel_function
	# We initalise the return value to be of the same type as x and to be constructed from 0.
	retval = type(x)(0)
	# Now we are going to iterate over the limit provided by the psl() method of x. Looking at the Bessel function
	# definition through power series, the starting degree will be n and the step will be 2.
	l = x.psl(n,2)
	for m in range(0,l):
		tmp = (x / 2)**(2 * m + n)
		tmp *= (-1)**m
		# Please note that the factorial function is defined in Pyranha but not in standard Python.
		tmp /= factorial(m) * factorial(m + n)
		retval += tmp
	return retval
