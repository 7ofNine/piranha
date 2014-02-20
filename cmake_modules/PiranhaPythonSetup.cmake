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

MACRO(PIRANHA_PYTHON_SETUP)
	# Find Python libraries
	FIND_PACKAGE(PythonLibs)
	IF(NOT PYTHONLIBS_FOUND)
		MESSAGE(FATAL_ERROR "No Python libraries found.")
	ENDIF(NOT PYTHONLIBS_FOUND)
	MESSAGE(STATUS "Python libraries: " "${PYTHON_LIBRARIES}")
	INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
	MESSAGE(STATUS "Python library: " "${PYTHON_LIBRARY}")
	IF(UNIX)
		STRING(REGEX MATCH libpython[0-9]\\.?[0-9] PYTHON_LIBRARY_VERSION_DOT ${PYTHON_LIBRARY})
		STRING(REGEX REPLACE libpython "" PYTHON_LIBRARY_VERSION_DOT ${PYTHON_LIBRARY_VERSION_DOT})
		STRING(REGEX REPLACE \\. "" PYTHON_LIBRARY_VERSION ${PYTHON_LIBRARY_VERSION_DOT})
		# In certain systems the Python lib is in /usr/lib(64), in others under /usr/lib/python2.5/config/.
		# We try here to catch both cases.
		STRING(REGEX REPLACE "(python${PYTHON_LIBRARY_VERSION_DOT}/config/)?libpython.*"
			"python${PYTHON_LIBRARY_VERSION_DOT}/site-packages/" PYTHON_MODULES_PATH ${PYTHON_LIBRARY})
	ELSE(UNIX)
		STRING(REGEX MATCH python[0-9]\\.?[0-9] PYTHON_LIBRARY_VERSION ${PYTHON_LIBRARY})
		STRING(REGEX REPLACE python "" PYTHON_LIBRARY_VERSION ${PYTHON_LIBRARY_VERSION})
#		SET(PYTHON_LIBRARY_VERSION ${PYTHONLIBS_VERSION_STRING})
#		SET(PYTHON_LIBRARY_VERSION "${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}")
		SET(PYTHON_MODULES_PATH .)
	ENDIF(UNIX)
	MESSAGE(STATUS "Python library version: " ${PYTHON_LIBRARY_VERSION})
	SET(PYDEXTENSION FALSE)
	IF(${PYTHON_LIBRARY_VERSION} GREATER 24 AND WIN32)
		MESSAGE(STATUS "Python >= 2.5 detected on WIN32 platform. Output extension for compiled modules will be '.pyd'.")
		SET(PYDEXTENSION TRUE)
	ENDIF(${PYTHON_LIBRARY_VERSION} GREATER 24 AND WIN32)
	MESSAGE(STATUS "Python modules install path: " "${PYTHON_MODULES_PATH}")
	FIND_PACKAGE(PythonInterp)
	IF(NOT PYTHONINTERP_FOUND)
		MESSAGE(FATAL_ERROR "No Python interpreter found.")
	ENDIF(NOT PYTHONINTERP_FOUND)
#	MESSAGE(STATUS "Python Interpretr: ${PYTHON_EXECUTABLE}")
	STRING(REPLACE "python.exe" "" PIRANHA_PYTHON_PATH ${PYTHON_EXECUTABLE})
	MESSAGE(STATUS "Python Path: ${PIRANHA_PYTHON_PATH}")
ENDMACRO(PIRANHA_PYTHON_SETUP)
