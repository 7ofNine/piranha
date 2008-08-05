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

# Try to find PyQt4 utilities, pyuic4 and pyrcc4:
# PYUIC4BINARY - Location of pyuic4 executable
# PYRCC4BINARY - Location of pyrcc4 executable
# PyQt4_FOUND - PyQt4 utilities found.

# Also provides macro similar to FindQt4.cmake's WRAP_UI and WRAP_RC,
# for the automatic generation of Python code from Qt4's user interface
# ('.ui') and resource ('.qrc') files. These macros are called:
# - PYQT4_WRAP_UI
# - PYQT4_WRAP_RC

IF(PYUIC4BINARY AND PYRCC4BINARY)
	# Already in cache, be silent
	SET(PyQt4_FIND_QUIETLY TRUE)
ENDIF(PYUIC4BINARY AND PYRCC4BINARY)

FIND_PROGRAM(PYUIC4BINARY pyuic4)
FIND_PROGRAM(PYRCC4BINARY pyrcc4)

MACRO(PYQT4_WRAP_UI outfiles)
	FOREACH(it ${ARGN})
		GET_FILENAME_COMPONENT(outfile ${it} NAME_WE)
		GET_FILENAME_COMPONENT(infile ${it} ABSOLUTE)
		SET(outfile ${CMAKE_CURRENT_BINARY_DIR}/ui_${outfile}.py)
		ADD_CUSTOM_TARGET(${it} ALL
			DEPENDS ${outfile}
		)
		ADD_CUSTOM_COMMAND(OUTPUT ${outfile}
			COMMAND ${PYUIC4BINARY} ${infile} -o ${outfile}
			MAIN_DEPENDENCY ${infile}
		)
		SET(${outfiles} ${${outfiles}} ${outfile})
	ENDFOREACH(it)
ENDMACRO (PYQT4_WRAP_UI)

MACRO(PYQT4_WRAP_RC outfiles)
	FOREACH(it ${ARGN})
		GET_FILENAME_COMPONENT(outfile ${it} NAME_WE)
		GET_FILENAME_COMPONENT(infile ${it} ABSOLUTE)
		SET(outfile ${CMAKE_CURRENT_BINARY_DIR}/${outfile}_rc.py)
		ADD_CUSTOM_TARGET(${it} ALL
			DEPENDS ${outfile}
		)
		ADD_CUSTOM_COMMAND(OUTPUT ${outfile}
			COMMAND ${PYRCC4BINARY} ${infile} -o ${outfile}
			MAIN_DEPENDENCY ${infile}
		)
		SET(${outfiles} ${${outfiles}} ${outfile})
	ENDFOREACH(it)
ENDMACRO (PYQT4_WRAP_RC)

IF(EXISTS ${PYUIC4BINARY} AND EXISTS ${PYRCC4BINARY})
   SET(PyQt4_FOUND TRUE)
ENDIF(EXISTS ${PYUIC4BINARY} AND EXISTS ${PYRCC4BINARY})

IF(PyQt4_FOUND)
	IF(NOT PyQt4_FIND_QUIETLY)
		MESSAGE(STATUS "Found PyQt4: ${PYUIC4BINARY}, ${PYRCC4BINARY}")
	ENDIF(NOT PyQt4_FIND_QUIETLY)
ELSE(PyQt4_FOUND)
	IF(PyQt4_FIND_REQUIRED)
		MESSAGE(FATAL_ERROR "Could NOT find PyQt4")
	ENDIF(PyQt4_FIND_REQUIRED)
ENDIF(PyQt4_FOUND)
