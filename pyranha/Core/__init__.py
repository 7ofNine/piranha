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

#print("Core.__init__.py.20 " + str(dir())) #DEBUG
from ._Core import *
from .impl import *

def is_iteratable(arg):
	"""
	Returns True if arg is iteratable, false otherwise.
	"""
	try:
		iter(arg)
		return True
	except TypeError:
		return False

def copy(arg):
	"""
	Standard copy function. Lifted from the copy module.
	"""
	import copy as __copy
	return __copy.copy(arg)

def latex(arg):
	"""
	Return latex-formatted string representation of the object.
	"""
	if isinstance(arg,int):
		return str(arg)
	if isinstance(arg,float):
		# Handle the case in which float representation is in scientific format.
		tmp = str(arg).split('e')
		assert(len(tmp) >= 1)
		if len(tmp) == 1:
			return tmp[0]
		retval = r'%s \cdot 10^{%s}' % (tmp[0],str(int(tmp[1])))
		return retval
	if isinstance(arg,complex):
		retval = r'\left( %s' % latex(arg.real)
		if arg.imag >= 0:
			retval += r'+'
		retval += r'%s\,\imath\right)' % latex(arg.imag)
		return retval
	if not hasattr(arg,"_latex_"):
		raise AttributeError('Object does not provide a _latex_() method.')
	return arg._latex_()

def latex_tab(series, width = .8 , geometry = 'a4paper,margin=0.2in', textsize = 'tiny'):
	"""
	Returns a table-like latex representation for an input series s. The geometry of the page (in latex geometry
	package format) and the latex textsize value can be passed as optional parameters. The width parameter, which must
	be in the interval ]0,1[, is the fraction of total space occupied by the first column, displaying the series' coefficients.
	The second column displays the series' keys.

	Argument s is expected to be a collection of indexable pairs, so that it is possible to pass not only a series but also, e.g., a list of
	coefficient-key pairs.

	In order to use the output of this function, the 'breqn', 'xtab', 'nicefrac' and 'geometry' packages must be present (they should
	be fairly common on most latex installations).

	PLEASE NOTE: apparently the xtab package may require many runs of the latex command in order to figure out exactly the placement of the long
	table resulting from a big series.
	"""
	if width <= 0 or width >= 1:
		raise ValueError('Use a width value in the ]0,1[ range.')
	l_w = str(float(width) * .9)
	r_w = str((1. - width) * .9)
	retval = r"""
	\documentclass[landscape]{article}
	\usepackage[""" + geometry + r"""]{geometry}
	\usepackage{xtab}
	\usepackage{breqn}
	\usepackage{nicefrac}
	\begin{document}\%s""" % textsize + r"""
	\begin{center}
	\begin{xtabular}{|p{""" + l_w + r"""\textwidth}|p{""" + r_w + r"""\textwidth}|} \hline
	"""
	for t in series:
		retval += r'\begin{minipage}{'+ l_w + r'\textwidth} \begin{dmath*} ' + latex(t[0]) + r' \end{dmath*} \end{minipage} & \begin{minipage}{' + r_w + r'\textwidth} \begin{dmath*}' \
			+ latex(t[1]) + r'\end{dmath*} \end{minipage} \\ \hline '
	retval += r"""
	\end{xtabular}
	\end{center}
	\end{document}

	"""
	return retval

def display(series):
	import shutil, tempfile, os.path
	from pyranha.Core import latex_tab
	from subprocess import Popen, PIPE, STDOUT
	if isinstance(series,str):
		s = series
	else:
		s = latex_tab(series.split())
	tmp_dir = tempfile.mkdtemp()
	try:
		tex_file = tempfile.NamedTemporaryFile(dir = tmp_dir, delete = False)
		tex_file.write(s)
		tex_file.close()
		proc = Popen(['pdflatex','-interaction=batchmode','-jobname=pyranha_output',tex_file.name], cwd = tmp_dir, stdout = PIPE, stderr = STDOUT)
		output = proc.communicate()[0]
		if proc.returncode:
			raise RuntimeError(output)
		while True:
			proc = Popen(['pdflatex','-interaction=batchmode','-jobname=pyranha_output',tex_file.name], cwd = tmp_dir, stdout = PIPE, stderr = STDOUT)
			tmp = proc.communicate()[0]
			if proc.returncode:
				raise RuntimeError(tmp)
			if tmp == output:
				break
			output = tmp
		proc = Popen(['okular',os.path.join(tmp_dir,'pyranha_output.pdf')], stdout = PIPE, stderr = STDOUT)
		proc.wait()
		input('Press Enter to continue...')
	finally:
		shutil.rmtree(tmp_dir)

psyms = impl.__psyms()
series = impl.__series()

def load(*args):
	"""
	Load series from list of filenames arguments.
	Type of series will be inferred - if possible - from the files' extensions.
	"""
	try:
		import IPython.ipapi
	except ImportError:
		raise ImportError("IPython not available.")
	ip = IPython.ipapi.get()
	args_list = list(args)
	for i in args_list:
		s = i.split(".")
		var_name = s.pop(0)
		if not s:
			raise ValueError("Could not establish an extension for the filename \"" + arg + "\".")
			return
		ext_name = s.pop()
		try:
			ip.ex("%s = %s(\"%s\")" % (var_name,ext_name,i))
		except NameError:
			raise NameError("Series type \"" + ext_name + "\" is unknown.")

settings = _Core.__settings()
stats = _Core.__stats()

def gui():
	try:
		import pyranha.Gui
		pyranha.Gui.mw.show()
	except ImportError:
		print("Gui support is not available or PyQt4 is not installed.")

#print("Core.__init__.py.173 " + str(dir()))   #DEBUG
