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

def alt_besselJ(o,x):
    """
    Python implementation of the power series expansion for Bessel functions of the first kind
    of order o and argument x.

    This implementation will give the same output as the besselJ() method of Pyranha's series
    (which is implemented in C++).
    """
    # Import the Pyranha module.
    import pyranha
    # First we check that o is an integer.
    if type(o) != int:
        raise TypeError("o must be integer.")
    # Check also that the type of x belongs to the available Pyranha manipulators.
    if type(x) not in pyranha.manipulators_type_tuple:
        raise TypeError("x is not a Pyranha series type.")
    # Next we check that o is positive. If not, change it to its negative.
    n = o
    if o < 0:
        n = -o
    # Now let's proceed to the actual expansion. Formula is available on Wikipedia, for reference:
    # http://en.wikipedia.org/wiki/Bessel_function
    # We initalise the return value to be of the same type as x and to be constructed from 0.
    retval = type(x)(0)
    # Now we are going to iterate over the limit provided by the psi() method of x. Looking at the Bessel function
    # definition through power series, the starting degree will be n and the step will be 2.
    l = x.psi(n,2)
    for m in range(0,l):
        tmp = (x / 2)**(2 * m + n)
        tmp *= (-1)**m
        # Here we use the unbound factorial functions of the series' type to avoid
        # losing precision when working with multi-precision coefficients.
        tmp *= (type(x).factorial(m) * type(x).factorial(m + n))**-1
        retval += tmp
    # If we corrected the order's sign we must also multiply by (-1)**n.
    if o < 0:
        retval *= (-1)**n
    return retval
