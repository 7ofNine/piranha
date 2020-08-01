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
# TODO: remove as many as possible of the exec function calls. BAd style and they have to be modified anyways
# because of different scoping effects from Python2 as a statemet (!!) vs. Python3 as function

# remove all the temporarily added print statements for debugging

# Most of your config files and extensions will probably start with this import
import os;
os.add_dll_directory("D:\\temp for tests\Piranhatest") # needed since 3.8 to be able to load dependency dll's

from IPython import get_ipython
ip = get_ipython()

#import IPython.ipapi
#ip = IPython.ipapi.get()

def main():
#   o = ip.options
#   o.system_verbose = 0
    #print("main.31")  #DEBUG 
    import pyranha    # this is where most of the import is done!
    #print("main.32")  #DEBUG
    for i in pyranha.__manipulators__:
        #print("main.35 " + i) #DEBUG 
        ip.ex("from pyranha import %s" % i)
        ip.ex("from pyranha.%s import %s" % (i,i.lower()))
        #print("main.38") #DEBUG
        # Try importing the complex counterpart.
        try:
            #print("main.41") #DEBUG 
            ip.ex("from pyranha.%s import %s" % (i,i.lower()+'c'))
        except:
            #print("main.44") #DEBUG
            pass
    for i in [x for x in pyranha.__all__ if x not in pyranha.__manipulators__]:
        #print("main.47: " + i) #DEBUG
        if i != "Gui" and i != "Test" and i != "Truncators":
            ip.ex("from pyranha.%s import *" % i)
    # Import default series type.
    #print("main.51") #DEBUG
    #print(dir()) #DEBUG
    #print(dir(pyranha)) #DEBUG
    ip.ex("from pyranha import ds")
    # Import test module.
    #print("main.54") #DEBUG
    #print(dir())
    ip.ex("from pyranha import Test")
    # Import truncator-handling class.
    #print("main.57") #DEBUG
    ip.ex("from pyranha import truncators")
    import_error_msg = """
        Warning: some of Pyranha's capabilities rely on numpy and matplotlib.
        Please consider installing these packages:
        http://numpy.scipy.org
        http://matplotlib.sf.net"""
    error_msg = False
    try:
        #print("main.66") #DEBUG
        ip.ex("import numpy")   # isn't the other way arround simpler? matplotlib relies on numpy?
        print("Numpy was successfully loaded.")
    except ImportError:
        if not error_msg: print(import_error_msg)
        error_msg = True
    try:
        #print("main.73") #DEBUG
        ip.ex("import matplotlib")
        ip.ex("matplotlib.interactive(True)")
        print("Matplotlib was successfully loaded. Interactive mode has been activated.")
    except ImportError:
        if not error_msg: print(import_error_msg)
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

if(__name__ == "__main__"):
    main()
