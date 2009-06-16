# -*- coding: iso-8859-1 -*-
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
		if i != "Gui" and i != "Test":
			ip.ex("from pyranha.%s import *" % i)
	# Import default series type.
	ip.ex("from pyranha import ds")
	# Import test module.
	ip.ex("from pyranha import Test")
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

def piranha_editor(self, filename, linenum=None):
	import os
	# Try to fetch the "EDITOR" environment variable, otherwise use the provided npp editor.
	try:
		editor = os.environ["EDITOR"]
	except KeyError:
		editor = 'npp\\notepad++.exe'
	# Marker for at which line to open the file (for existing objects)
	if linenum is None or editor=='notepad' or editor=='npp\\notepad++.exe':
		linemark = ''
	else:
		linemark = '+%d' % int(linenum)
	# Enclose in quotes if necessary and legal
	if ' ' in editor and os.path.isfile(editor) and editor[0] != '"':
		editor = '"%s"' % editor
	# Call the actual editor
	os.system('%s %s %s' % (editor,linemark,filename))

ip.set_hook('editor', piranha_editor)

main()
