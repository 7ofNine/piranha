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

from @MODULE_NAME@ import *
import re

try:
  from pylab import plot,semilogy
  def __plot_range_evaluator(self,*args,**kwargs):
    plot(self.times(),self.values(),*args,**kwargs)
  def __plot_tc(self,*args,**kwargs):
    values = list()
    for i in range(self.plain_eval().times().__len__()):
      values.append(abs(self.plain_eval().values()[i]-self.manip_eval().values()[i]))
    semilogy(self.plain_eval().times(),self.plain_eval().values(),*args,**kwargs)
    semilogy(self.manip_eval().times(),values,*args,**kwargs)
  range_evaluator.plot = __plot_range_evaluator
  tc_list=[i for i in dir(@MODULE_NAME@) if re.search('^tc_.*',i)]
  for i in tc_list:
    eval(i).plot = __plot_tc
except ImportError:
  print "Matplotlib is not installed, disabling plotting methods."
except NameError:
  pass
