""" User configuration file for IPython

This is a more flexible and safe way to configure ipython than *rc files
(ipythonrc, ipythonrc-pysh etc.)

This file is always imported on ipython startup. You can import the
ipython extensions you need here (see IPython/Extensions directory).

Feel free to edit this file to customize your ipython experience.

Note that as such this file does nothing, for backwards compatibility.
Consult e.g. file 'ipy_profile_sh.py' for an example of the things 
you can do here.

See http://ipython.scipy.org/moin/IpythonExtensionApi for detailed
description on what you could do here.
"""

# Most of your config files and extensions will probably start with this import

import IPython.ipapi
ip = IPython.ipapi.get()

def main():
	o = ip.options
	o.system_verbose = 0
	ip.ex("import math")
	ip.ex("import pyranha")
	import pyranha
	for i in pyranha.__manipulators__:
		ip.ex("from pyranha import %s" % i)
		ip.ex("from pyranha.%s import %s" % (i,i.lower()))
		# Try importing the complex counterpart.
		try:
			ip.ex("from pyranha.%s import %s" % (i,i.lower()+'c'))
		except:
			pass
	for i in filter(lambda x: x not in pyranha.__manipulators__,pyranha.__all__):
		if i != "Gui":
			ip.ex("from pyranha.%s import *" % i)
	ip.ex("from pyranha import ds")
	import_error_msg = """
		Warning: many of Pyranha's capabilities rely on numpy and matplotlib.
		Please consider installing these packages:
		http://numpy.scipy.org
		http://matplotlib.sf.net"""
	error_msg = False
	try:
		ip.ex("import numpy")
		print "Numpy was successfully loaded."
	except ImportError:
		if not error_msg: print import_error_msg
		error_msg = True
	try:
		ip.ex("import matplotlib")
		ip.ex("matplotlib.interactive(True)")
		print "Matplotlib was successfully loaded. Interactive mode has been activated."
	except ImportError:
		if not error_msg: print import_error_msg
		error_msg = True
	ip.ex("settings.format = pretty")

main()
