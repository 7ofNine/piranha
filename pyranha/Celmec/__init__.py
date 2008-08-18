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
