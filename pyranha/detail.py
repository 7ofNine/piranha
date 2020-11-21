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


# TODO: get rid of exec !

def __build_manipulators(manipulators):
    """
    Build the list of manipulator types
    """
    #print("pyranha.detail.py.24")     #DEBUG
    type_list = []
    for i in manipulators:
        #print("pyranha.detail.py.30: " + str(i))   #DEBUG
        __cur_manip = {} # just a work arround before we can get rid of the exec function calls
        exec("import pyranha.%s as __cur_manip" % i , {}, __cur_manip)
        cur_manip = __cur_manip['__cur_manip']  # get the actual object
        type_list.append(getattr(cur_manip, i.lower()))
        # Let's try to see if we can get the complex counterpart.
        try:
            type_list.append(getattr(cur_manip,i.lower()+'c'))
        except AttributeError:
            pass
    #print("pyranha.detail.py.34 " + str(type_list))   #DEBUG
    return tuple(type_list)

def __series_filtered(self, criterion = None):
    """
    Return input series with terms filtered according to criterion.

    If criterion is a callable binary object, the series will be decomposed into a list of coefficient-key pairs and the callable
    will be invoked using coefficients and keys as arguments. If the return value of the callable is True, then
    coefficient and key are retained, otherwise they are eliminated from the series.

    If criterion is a list of callable binary objects, each callable will be called recursively to filter out
    coefficient-key pairs, beginning with the current series and descending in a recursive fashion into the coefficient
    series (and from there into the coefficient series of the coefficient series, etc.).

    A None criterion is interpreted as a callable that always returns True.

    Please note that since reassembling coefficient-key pairs involves series multiplications, active truncators
    have an effect on the filtering process.
    """
    #print("pyranha.detail.py.52")   #DEBUG
    from copy import copy
    if criterion is None:
        return copy(self)
    try:
        iter(criterion)
        crit = list(criterion)
    except TypeError:
        crit = [criterion]
    # If list is empty, return self.
    if not crit:
        return copy(self)
    retval = type(self)()
    rec_depth = len(crit)
    if rec_depth > len(self.arguments):
        raise ValueError('Cannot apply %d recursive filters to a series of echelon level %d' % (rec_depth,len(self.arguments) - 1))
    for c in crit:
        if not callable(c) and not c is None:
            raise ValueError('Please provide one or more unary callables (or None) as filtering criterions.')
    def filter_series(cur,tot,s,crits):
        if cur > tot:
            return copy(s)
        retval = type(s)()
        l = s.split(cur)
        for t in l:
            if crits[cur] is None or crits[cur](t):
                retval += filter_series(cur + 1, tot, t[0], crits) * t[1]
        return retval
    return filter_series(0,rec_depth - 1,self,crit)

def __series_short_type(self):
    """
    Return a short string containing the series' class name.
    """
    #print("pyranha.detail.py.85")  #DEBUG
    return str(type(self)).rpartition('.')[-1].strip('>\'')

def __series_psi(self, n = 0, s = 1):
    """
    Return the limit of the power series expansion compatible with the truncator currently in use by the series.
    The optional arguments are the starting degree (n) and the step size (s) of the power series expansion.
    """
    #print("pyranha.detail.py.92")  #DEBUG
    return self.__psi__(n,s)

def __series_split(self, n = 0):
    """
    Split the series into a sequence of coefficient-key pairs.

    For n > 0, this method will try to split the series at higher echelon levels. That is,
    if n == 1 and the series is degenerate (i.e., one single term with unitary key), a split()
    on the single term's coefficient series will take place and the return value will be a
    sequence of coefficient-key pairs for the coefficient series.
    """
    #print("pyranha.detail.py.103") #DEBUG
    return self.__split__(n)

def __series_eval(self,arg):
    """
    Numerical evaluation of the series.

    If arg is a floating point number, then arg is assumed to be a time and the series is evaluated according to the
    time evaluation vectors of each psym of the series.

    If arg is a dictionary of string-floating point value pairs, the series is evaluated as if each psyms were substituted the numerical
    value specified in the dictionary. If the dictionary does not contain an entry for each psym of the series, a ValueError exception
    will be raised.

    If arg is anything else, a TypeError exception will be raised.
    """
    #print("pyranha.detail.py.118") #DEBUG
    if isinstance(arg,(float,int)):
        return self.__eval__(arg)
    if isinstance(arg,dict):
        from pyranha.Core import EvalDict
        d = EvalDict()
        for t in arg:
            # Check that all entries are of the right type.
            if not isinstance(t,str) or not isinstance(arg[t],(float,int)):
                raise TypeError('The provided dictionary does not contain only string-float pairs.')
            d[t] = arg[t]
        return self.__eval__(d)
    raise TypeError('Cannot use the type ' + str(type(arg)) + ' for evaluation.')

def __series_repr(self):
    """
    __repr__ method that prints the series' type and then the series itself in pretty print.
    """
    #print("pyranha.detail.py.135") #DEBUG
    retval = '%s series: ' % self.__short_type__
    retval += self.__impl_repr__()
    return retval

def __series_contains(self,name):
    """
    Check whether symbol named 'name' is present in the series.
    """
    #print("pyranha.detail.py.143") #DEBUG
    from pyranha.Core import psym
    return psym(name) in reduce(lambda a,b: a + b, self.arguments)

def __series_deepcopy(self,memo):
    #print("pyranha.detail.py.147") #DEBUG
    return self.__copy__()

def __add_method(module_name, method_name, function):
    """
    Add a method to a manipulator.
    """
    loc = {'module_name' : module_name, 'method_name' : method_name, 'function' : function} 
    #print("pyranha.detail.py.153 " + module_name ) #DEBUG
    exec("import %s.%s as __cur_manip" % ("pyranha", module_name), {}, loc)
    #print("pyranha.detail.py.155 " + str(__cur_manip)) #DEBUG
    exec("__cur_manip.%s.%s = function" % (module_name.lower(), method_name), {}, loc)
    # Try to take care of the complex counterpart.
    try:
        exec("__cur_manip.%s.%s = function " % ((module_name.lower()+'c'), method_name), {}, loc)
    except AttributeError:
        pass

def __add_property(module_name, property_name, fget=None, fset=None, fdel=None, doc=None):
    """
    Add a property to a manipulator.
    """
    loc = {'module_name' : module_name, 'property_name' : property_name, 'fget' : fget, 'fset' : fset, 'fdel' : fdel, 'doc' : doc }
    #print("pyranha.detail.py.165 " + module_name) #DEBUG
    exec("import %s.%s as __cur_manip" % ("pyranha", module_name), {}, loc )   # how can we avoid exec !!! causes issues
    #print("pyranha.detail.py.167 :" + str(__cur_manip['__cur_manip']) + '\n') #DEBUG #exec keeps __cur_manip it locally and it is not visible. That is Python3!!
    exec("__cur_manip.%s.%s = property(fget,fset,fdel,doc)" % (module_name.lower(), property_name), {}, loc)
    try:
        #print("pyranha.detail.py.187")  #DEBUG
        exec("__cur_manip.%s.%s = property(fget,fset,fdel,doc) " % ((module_name.lower()+'c'),property_name), {}, loc)
    except AttributeError:
        #print("pyranha.detail.py.190")  #DEBUG
        pass

def __enhance_manipulators(manipulators):
    """
    Add useful methods and properties to a list of manipulators.
    """
    for i in manipulators:
        #print("pyranha.detail.py.177 " + str(dir(i)))  #DEBUG
        __add_property(i, "__short_type__", __series_short_type)
        __add_method(i, "filtered", __series_filtered)
        __add_method(i, "psi", __series_psi)
        __add_method(i, "split", __series_split)
        __add_method(i, "eval", __series_eval)
        __add_method(i, "__contains__", __series_contains)
        __add_method(i, "__repr__", __series_repr)
        __add_method(i, "__deepcopy__", __series_deepcopy)
        #print("pyranha.detail.py.186 " + str(dir(i)))  #DEBUG

#print("pyranha.detail.py.186a")
import pyranha   # if not imported it doesn't know the name pyranha for pyranha.__all__. Why???doesn't seem to do anything. It's a recursive call. This already called during an import pyranha!!
#print("pyranha.detail.py.188")
from functools import reduce

# Let's try do define a default series type
#print(__name__)  # for test only scope is pyranha.detail
#print("Trace 2")
if len([x for x in pyranha.__all__ if x not in pyranha.__manipulators__]) > 0:
    #print("pyranha.detail.py.192")
    #exec ("import %s as __last_manipulator" % pyranha.__manipulators__[-1])
    loc = {}
    glob = {'__name__' : __name__}
    exec("from . import %s as __last_manipulator" % pyranha.__manipulators__[-1], glob, loc)
    last_manipulator = loc['__last_manipulator']
    #print("pyranha.detail.py.194" + str(last_manipulator)) #DEBUG
    ds = getattr(last_manipulator, pyranha.__manipulators__[-1].lower(), "") # fix it: for setting the default series type. worked in Python 2
    #print("pyranha.detail.py.195") #DEBUG
    setattr(pyranha, 'ds', ds)
    #print("pyranha.detail.py.198") #DEBUG
    print(("Default series type is " + str(ds)))
else:
    print("Default series type could not be established, assigning None.")
    ds = None

#print("pyranha.detail.py.204") #DEBUG  
__enhance_manipulators(pyranha.__manipulators__)

#print("pyranha.detail.py.207") #DEBUG
manipulators = tuple(__build_manipulators(pyranha.__manipulators__))
#print("pyranha.detail.py.236") #DEBUG