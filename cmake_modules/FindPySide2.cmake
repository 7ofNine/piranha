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

# Try to find PySide2 utilities, pyuic4 and pyrcc4:
# PYSIDE2UICBINARY - Location of pyuic4 executable
# PYSIDE2RCBINARY - Location of pyrcc4 executable
# PySIDE2_FOUND - pySide2 utilities found.

# Also provides macro similar to FindQt4.cmake's WRAP_UI and WRAP_RC,
# for the automatic generation of Python code from Qt4's user interface
# ('.ui') and resource ('.qrc') files. These macros are called:
# - PYSIDE2_WRAP_UI
# - PYSIDE2_WRAP_RC

IF(PYSIDE2UICBINARY AND PYSIDE2RCCBINARY)
    # Already in cache, be silent
    SET(PySide2_FIND_QUIETLY TRUE)
ENDIF(PYSIDE2UICBINARY AND PYSIDE2RCCBINARY)

FIND_PROGRAM(PYSIDE2UICBINARY pyside2-uic)
FIND_PROGRAM(PYSIDE2RCCBINARY pyside2-rcc)

MACRO(PYSIDE2_WRAP_UI outfiles)
    FOREACH(it ${ARGN})
        GET_FILENAME_COMPONENT(outfile ${it} NAME_WE)
        GET_FILENAME_COMPONENT(infile ${it} ABSOLUTE)
        SET(outfile ${CMAKE_CURRENT_BINARY_DIR}/ui_${outfile}.py)
        ADD_CUSTOM_TARGET(${it} ALL
            DEPENDS ${outfile}
        )
        ADD_CUSTOM_COMMAND(OUTPUT ${outfile}
            COMMAND ${PYSIDE2UICBINARY} ${infile} -o ${outfile}
            MAIN_DEPENDENCY ${infile}
        )
        SET(${outfiles} ${${outfiles}} ${outfile})
    ENDFOREACH(it)
ENDMACRO (PYSIDE2_WRAP_UI)

MACRO(PYSIDE2_WRAP_RC outfiles)
    FOREACH(it ${ARGN})
        GET_FILENAME_COMPONENT(outfile ${it} NAME_WE)
        GET_FILENAME_COMPONENT(infile ${it} ABSOLUTE)
        SET(outfile ${CMAKE_CURRENT_BINARY_DIR}/${outfile}_rc.py)
        ADD_CUSTOM_TARGET(${it} ALL
            DEPENDS ${outfile}
        )
        ADD_CUSTOM_COMMAND(OUTPUT ${outfile}
            COMMAND ${PYSIDE2RCCBINARY} ${infile} -o ${outfile}
            MAIN_DEPENDENCY ${infile}
        )
        SET(${outfiles} ${${outfiles}} ${outfile})
    ENDFOREACH(it)
ENDMACRO (PYSIDE2_WRAP_RC)

IF(EXISTS ${PYSIDE2UICBINARY} AND EXISTS ${PYSIDE2RCCBINARY})
   SET(PySide2_FOUND TRUE)
ENDIF(EXISTS ${PYSIDE2UICBINARY} AND EXISTS ${PYSIDE2RCCBINARY})

IF(PySide2_FOUND)
    IF(NOT PySide2_FIND_QUIETLY)
        MESSAGE(STATUS "Found PySide2: ${PYSIDE2UICBINARY}, ${PYSIDE2RCCBINARY}")
    ENDIF(NOT PySide2_FIND_QUIETLY)
ELSE(PySide2_FOUND)
    IF(PySide2_FIND_REQUIRED)
        MESSAGE(FATAL_ERROR "Could NOT find PySide2")
    ENDIF(PySide2_FIND_REQUIRED)
ENDIF(PySide2_FOUND)
