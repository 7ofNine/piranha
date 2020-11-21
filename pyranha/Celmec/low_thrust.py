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

__doc__ = 'Theories of keplerian motion perturbed by constant thrust in a fixed direction.'

import pyranha.Celmec

class lt_base(object):
    def __init__(self,thrust,axis,trunc_order):
        if not isinstance(thrust,float) or thrust <= 0:
            raise ValueError('Invalid thrust.')
        if not isinstance(axis,int) or not axis in [0,1,2]:
            raise ValueError('Invalid axis.')
        if not isinstance(trunc_order,int) or trunc_order < 3:
            raise ValueError('Invalid truncation order.')
        self.__thrust = thrust
        self.__axis = axis
        self.__trunc_order = trunc_order

    def _numerical_integrator(self,s0,tf,n):
        """
        Return a list of n state vectors starting from t = 0 to t = tf assuming as initial conditions
        the state s0 (list of 6 numbers, rectangular position + velocity).
        """
        from scipy.integrate import odeint
        from numpy import linspace
        from numpy.linalg import norm
        axis = self.axis
        thrust = self.thrust
        def dyn(s,_):
            r = norm(s[:3])
            r3 = r * r * r
            retval = [s[3],s[4],s[5],-s[0] / r3,-s[1] / r3,-s[2] / r3]
            retval[3 + axis] = retval[3 + axis] + thrust
            return retval
        return odeint(dyn,s0,linspace(0,tf,n))

    def _H1(self,mdelaunay):
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

    def Hamiltonian(self,cartesian):
        from math import sqrt
        x, y, z, vx, vy, vz = cartesian
        v2 = vx ** 2 + vy ** 2 + vz ** 2
        r = sqrt(x ** 2 + y ** 2 + z ** 2)
        return v2 / 2. - 1. / r - self.thrust * cartesian[self.axis]

    @property
    def thrust(self):
        return self.__thrust

    @property
    def axis(self):
        return self.__axis

    @property
    def trunc_order(self):
        return self.__trunc_order

class lt_divergent_mdel(pyranha.Celmec.lie_theory,lt_base):
    def __init__(self,order,t,thrust,axis = 0, trunc_order = 5, verbose = True):
        from pyranha.Core import psym, rational
        from pyranha import manipulators, truncators
        if not t in manipulators:
            raise TypeError('t must be a series type.')
        if not isinstance(order,int) or order <= 0:
            raise ValueError('Invalid order.')
        # First constructor.
        lt_base.__init__(self,thrust,axis,trunc_order)
        Lam, P, Q, lam, p, q = t(psym('Lam')), t(psym('P')), t(psym('Q')), t(psym('lam')), t(psym('p')), t(psym('q'))
        lam0, dlam = t(psym('lam0')), t(psym('dlam'))
        eps_series = t(psym('eps'))
        # Setup the truncators.
        truncators.unset()
        truncators.degree.set(['P','Q'],rational(self.trunc_order,2))
        H1 = self._H1([Lam,P,Q,lam,p,q])
        truncators.unset()
        H = - (2 * Lam ** 2) ** -1 + eps_series * H1.sub('lam',lam0 + dlam)
        truncators.degree.set(['P','Q'],rational(self.trunc_order,2))
        # Second constructor.
        pyranha.Celmec.lie_theory.__init__(self,H,'eps',['Lam','P','Q'],['dlam','p','q'],[lambda Hn: Lam ** 3 * Hn.integrate('dlam')] * order,verbose = verbose)
        # Reset the truncator on exit.
        truncators.unset()

    def solve_last(self,init,t):
        from copy import deepcopy
        retval = deepcopy(init)
        retval['dlam'] = retval['dlam'] + 1. / (retval['Lam'] ** 3) * t
        return retval

    def numerical(self,tf,n):
        from pyranha.Celmec import mdelaunay2oe, oe2s
        init_list = self.init_list[0]
        s0 = oe2s(mdelaunay2oe([init_list['Lam'],init_list['P'],init_list['Q'],init_list['dlam'] + init_list['lam0'],init_list['p'],init_list['q']]))
        return self._numerical_integrator(s0,tf,n)

    def cartesian(self,t):
        from pyranha.Celmec import mdelaunay2oe, oe2s
        tmp = self.evaluate(t)
        return oe2s(mdelaunay2oe([tmp['Lam'],tmp['P'],tmp['Q'],tmp['dlam'] + tmp['lam0'],tmp['p'],tmp['q']]))

class lt_convergent_planar(pyranha.Celmec.lie_theory,lt_base):
    def __init__(self,order,t,thrust,axis = 0, trunc_order = 5, verbose = True):
        # TODO: fix for p > pi.
        from pyranha.Core import psym, rational
        from pyranha import manipulators, truncators
        if not axis in [0,1]:
            raise ValueError('Invalid axis.')
        if not t in manipulators:
            raise TypeError('t must be a series type.')
        if not isinstance(order,int) or order <= 0:
            raise ValueError('Invalid order.')
        # First ctor.
        lt_base.__init__(self,thrust,axis,trunc_order)
        Lam, P, Q, lam, p, q = t(psym('Lam')), t(psym('P')), t(psym('Q')), t(psym('lam')), t(psym('p')), t(psym('q'))
        eps_series = t(psym('eps'))
        # Setup the truncators.
        truncators.unset()
        truncators.degree.set(['P','Q'],rational(self.trunc_order,2))
        H = - (2 * Lam ** 2) ** -1 + eps_series * self._H1([Lam,P,Q,lam,p,q]).sub('Q',t(0))
        # Second ctor.
        pyranha.Celmec.lie_theory.__init__(self,H,'eps',['Lam','P','Q'],['lam','p','q'],[lambda Hn: Lam ** 3 * sum([t for t in Hn if t.h_degree('lam') != 0]).integrate('lam')] * order, verbose = verbose)
        # Isolate c1 and its derivatives.
        self.__c_list = [sum([t for t in self.H[i + 1] if t.degree('eps') == i + 1]).sub('eps',self.series_type(1)).sub('p',self.series_type(0)) for i in range(self.order)]
        self.__c_Lam_list = [c.partial('Lam') for c in self.__c_list]
        self.__c_P_list = [c.partial('P') for c in self.__c_list]
        # Reset the truncator on exit.
        truncators.unset()

    def cosp_P(self,P,init):
        from copy import deepcopy
        if self.order != 1:
            raise ValueError('Non-unitary order is not allowed.')
        H = self.H[1].eval(init)
        H0 = -1. / (2. * init['Lam'] ** 2)
        init_copy = deepcopy(init)
        init_copy['P'] = P
        return (H - H0) / (self.thrust * self.__c_list[0].eval(init_copy))

    def solve_last(self,init,t):
        from copy import deepcopy
        from math import sqrt
        from scipy.integrate import quad
        from scipy.optimize import fsolve
        if self.order != 1:
            raise ValueError('Cannot solve for non unitary order.')
        # Value of the hamiltonian in the last set of variables.
        H_last = self.H[1].eval(init)

        # Define functions for evaluation of c1(P) and its derivatives.
        def c1_P(P):
            init_copy = deepcopy(init)
            init_copy['P'] = P
            return self.__c_list[0].eval(init_copy)

        def c1_Lam_P(P):
            init_copy = deepcopy(init)
            init_copy['P'] = P
            return self.__c_Lam_list[0].eval(init_copy)

        def c1_P_P(P):
            init_copy = deepcopy(init)
            init_copy['P'] = P
            return self.__c_P_list[0].eval(init_copy)

        # cos(p) as a function of P.
        def cosp_P(P):
            return (H_last + 1. / (2. * init['Lam'] ** 2)) / (self.thrust * c1_P(P))

        # Define integrand function.
        def I_P(P):
            c1 = c1_P(P)
            c2p = cosp_P(P) ** 2
            return 1. / (self.thrust * c1 * sqrt(1. - c2p))

        # Time as a function of P.
        def t_P(P):
            return quad(I_P,init['P'],P)[0]

        # P as function of time.
        def P_t(t):
            return fsolve(lambda P: t_P(P) - t,init['P'])

        # cos(p) as a function of time.
        def cosp_t(t):
            return cosp_P(P_t(t))

        def lam_t(t):
            def f(time):
                P = P_t(time)
                return c1_Lam_P(P) * cosp_P(P)
            return init['lam'] + 1. / (init['Lam'] ** 3) * t + self.thrust * quad(f,0,t)[0]

        def p_t(t):
            # TODO: no need to integrate here... acos with checks?
            def f(time):
                P = P_t(time)
                return c1_P_P(P) * cosp_P(P)
            return init['p'] + self.thrust * quad(f,0,t)[0]

        # Build return values.
        init_copy = deepcopy(init)
        init_copy['P'] = P_t(t)
        init_copy['Q'] = 0.
        init_copy['lam'] = lam_t(t)
        init_copy['p'] = p_t(t)
        init_copy['q'] = 0.
        return init_copy

    def numerical(self,tf,n):
        from pyranha.Celmec import mdelaunay2oe, oe2s
        init_list = self.init_list[0]
        s0 = oe2s(mdelaunay2oe([init_list['Lam'],init_list['P'],0,init_list['lam'],init_list['p'],0]))
        return self._numerical_integrator(s0,tf,n)

    def cartesian(self,t):
        from pyranha.Celmec import mdelaunay2oe, oe2s
        tmp = self.evaluate(t)
        return oe2s(mdelaunay2oe([tmp['Lam'],tmp['P'],tmp['Q'],tmp['lam'],tmp['p'],tmp['q']]))

class lt_divergent_poincare(pyranha.Celmec.lie_theory,lt_base):
    def __init__(self,order,t,thrust,axis = 0, trunc_order = 5, verbose = True):
        from pyranha.Core import psym, rational
        from pyranha import manipulators, truncators
        if not t in manipulators:
            raise TypeError('t must be a series type.')
        if not isinstance(order,int) or order <= 0:
            raise ValueError('Invalid order.')
        # First constructor.
        lt_base.__init__(self,thrust,axis,trunc_order)
        Lam, P, Q, lam, p, q = t(psym('Lam')), t(psym('P')), t(psym('Q')), t(psym('lam')), t(psym('p')), t(psym('q'))
        lam0, dlam = t(psym('lam0')), t(psym('dlam'))
        eps_series = t(psym('eps'))
        two, x, y, z, v = t(psym('two')), t(psym('x')), t(psym('y')), t(psym('z')), t(psym('v'))
        # Complex type.
        ct = type(t().complex())
        # Setup the truncators.
        truncators.unset()
        truncators.degree.set(['P','Q'],rational(self.trunc_order,2))
        H1 = self._H1([Lam,P,Q,lam,p,q])
        truncators.unset()
        H = - (2 * Lam ** 2) ** -1 + eps_series * H1.sub('lam',lam0 + dlam)
        # Exact conversion to Poincare' variables.
        H = H.ei_sub('p',ct(y * (two * P) ** (-rational(1,2)),x * (two * P) ** (-rational(1,2)))).ei_sub('q',ct(z * (two * Q) ** (-rational(1,2)),v * (two * Q) ** (-rational(1,2)))).sub('P',(x ** 2 + y ** 2) * two ** -1).sub('Q',(z ** 2 + v ** 2) * two ** -1)
        # Re-establish the truncator, this time on Poincare' variables.
        truncators.degree.set(['x','y','z','v'],rational(self.trunc_order,2))
        # Second constructor.
        pyranha.Celmec.lie_theory.__init__(self,H,'eps',['Lam','y','z'],['dlam','x','v'],[lambda Hn: Lam ** 3 * Hn.integrate('dlam')] * order, verbose = verbose)
        # Reset truncator on exit.
        truncators.unset()

    def solve_last(self,init,t):
        from copy import deepcopy
        retval = deepcopy(init)
        retval['dlam'] = retval['dlam'] + 1. / (retval['Lam'] ** 3) * t
        return retval

    def numerical(self,tf,n):
        from pyranha.Celmec import poincare2mdelaunay, mdelaunay2oe, oe2s
        init_list = self.init_list[0]
        s0 = oe2s(mdelaunay2oe(poincare2mdelaunay([init_list['Lam'],init_list['y'],init_list['z'],init_list['dlam'] + init_list['lam0'],init_list['x'],init_list['v']])))
        return self._numerical_integrator(s0,tf,n)

    def cartesian(self,t):
        from pyranha.Celmec import poincare2mdelaunay, mdelaunay2oe, oe2s
        tmp = self.evaluate(t)
        return oe2s(mdelaunay2oe(poincare2mdelaunay([tmp['Lam'],tmp['y'],tmp['z'],tmp['dlam'] + tmp['lam0'],tmp['x'],tmp['v']])))

class lt_convergent(pyranha.Celmec.lie_theory,lt_base):
    def __init__(self,order,t,thrust,trunc_order = 5,verbose = True):
        from pyranha.Core import psym, rational
        from pyranha import manipulators, truncators
        if not t in manipulators:
            raise TypeError('t must be a series type.')
        if not isinstance(order,int) or order <= 0:
            raise ValueError('Invalid order.')
        # First ctor.
        lt_base.__init__(self,thrust,2,trunc_order)
        Lam, P, Q, lam, p, q = t(psym('Lam')), t(psym('P')), t(psym('Q')), t(psym('lam')), t(psym('p')), t(psym('q'))
        eps_series = t(psym('eps'))
        # Setup the truncators.
        truncators.unset()
        truncators.degree.set(['P','Q'],rational(self.trunc_order,2))
        H = - (2 * Lam ** 2) ** -1 + eps_series * self._H1([Lam,P,Q,lam,p,q])
        # Second ctor.
        pyranha.Celmec.lie_theory.__init__(self,H,'eps',['Lam','P','Q'],['lam','p','q'],[lambda Hn: Lam ** 3 * sum([t for t in Hn if t.h_degree('lam') != 0]).integrate('lam')] * order,verbose = verbose)
        # Reset the truncator on exit.
        truncators.unset()

    def solve_last(self,init,t):
        raise NotImplementedError('Solution of this theory has not been implemented.')
