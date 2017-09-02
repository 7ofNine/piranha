# Copyright (C) 2007, 2008 by Francesco Biscani; 2017 Hartmuth Gutsche
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

# on MS Windows install Python with update to the registry in order for this to work. Pretty bad if you need to maintain 
# several versions

MACRO(PIRANHA_PYTHON_SETUP)
	FIND_PACKAGE(PythonInterp)
	IF(NOT PYTHONINTERP_FOUND)
		MESSAGE(FATAL_ERROR "No Python interpreter found.")
	ENDIF(NOT PYTHONINTERP_FOUND)
	STRING(REPLACE "python.exe" "" PIRANHA_PYTHON_PATH ${PYTHON_EXECUTABLE})
	MESSAGE(STATUS "Python Path:                 " ${PIRANHA_PYTHON_PATH})

	# Find Python libraries
	FIND_PACKAGE(PythonLibs)
	IF(NOT PYTHONLIBS_FOUND)
		MESSAGE(FATAL_ERROR "No Python libraries found.")
	ENDIF(NOT PYTHONLIBS_FOUND)
        
	MESSAGE(STATUS "Python libraries:            " ${PYTHON_LIBRARIES})
#	INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
        MESSAGE(STATUS "Python headers include path: " ${PYTHON_INCLUDE_PATH})
	MESSAGE(STATUS "Python library:              " ${PYTHON_LIBRARY})
        STRING(REGEX MATCH python[0-9]\\.?[0-9] PYTHON_LIBRARY_VERSION ${PYTHON_LIBRARY})
	STRING(REGEX REPLACE python "" PYTHON_LIBRARY_VERSION ${PYTHON_LIBRARY_VERSION})
        SET(PYTHON_MODULES_PATH .)
	MESSAGE(STATUS "Python library version:      " ${PYTHON_LIBRARY_VERSION})
	
	IF(${PYTHON_LIBRARY_VERSION} LESS 25 AND WIN32)
		MESSAGE(FATAL_ERROR STATUS "Python < 2.5 detected on WIN32 platform. This is not supported")
	ENDIF()
        
	MESSAGE(STATUS "Python modules install path: " "${PYTHON_MODULES_PATH}")
        
ENDMACRO(PIRANHA_PYTHON_SETUP)
