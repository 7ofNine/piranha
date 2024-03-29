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
    FIND_PACKAGE(Python3 COMPONENTS Interpreter Development ) 
    IF(NOT Python3_FOUND)
        MESSAGE(FATAL_ERROR "  -- No Python 3 interpreter found.--")
    else()
        SET(PIRANHA_PYTHON_PATH "${Python3_EXECUTABLE}")
        CMAKE_PATH(GET PIRANHA_PYTHON_PATH PARENT_PATH PIRANHA_PYTHON_ROOT)
        CMAKE_PATH(GET PIRANHA_PYTHON_ROOT FILENAME PIRANHA_PYTHON_BASE)
        STRING(TOLOWER ${PIRANHA_PYTHON_BASE} PIRANHA_PYTHON_BASE)
        SET(PIRANHA_PYTHON_LIBRARY_RELEASE "${Python3_LIBRARY_RELEASE}")
        SET(PIRANHA_PYTHON_LIBRARY_DEBUG "${Python3_LIBRARY_DEBUG}")
        SET(PIRANHA_PYTHON_LIBRARIES "${Python3_LIBRARIES}")
        SET(PIRANHA_PYTHON_LIBRARIES_DIR "${Python3_LIBRARY_DIRS}")
        SET(PIRANHA_PYTHON_INCLUDE_DIR "${Python3_INCLUDE_DIRS}")
        SET(PIRANHA_PYTHON_LIBRARY_VERSION ${Python3_VERSION})
        SET(PIRANHA_PYTHON_VERSION_MAJOR ${Python3_VERSION_MAJOR})
        CMAKE_PATH(APPEND PIRANHA_PYTHON_ROOT "Scripts" OUTPUT_VARIABLE PIRANHA_PYTHON_SCRIPTS)
    endif()

    MESSAGE(STATUS "Python Root:                 " ${PIRANHA_PYTHON_ROOT})
    MESSAGE(STATUS "Python Base:                 " ${PIRANHA_PYTHON_BASE})
	MESSAGE(STATUS "Python Path:                 " ${PIRANHA_PYTHON_PATH})
    MESSAGE(STATUS "Python headers include path: " ${PIRANHA_PYTHON_INCLUDE_DIR})
	MESSAGE(STATUS "Python library (Release)   : " ${PIRANHA_PYTHON_LIBRARY_RELEASE})
    MESSAGE(STATUS "Python library (Debug)     : " ${PIRANHA_PYTHON_LIBRARY_DEBUG})
    MESSAGE(STATUS "Python libraries           : " ${PIRANHA_PYTHON_LIBRARIES})
    MESSAGE(STATUS "Python libraries path      : " ${PIRANHA_PYTHON_LIBRARIES_DIR})
    MESSAGE(STATUS "Python library version     : " ${PIRANHA_PYTHON_LIBRARY_VERSION})
    MESSAGE(STATUS "Python scripts             : " ${PIRANHA_PYTHON_SCRIPTS})
    
	
    SET(PYTHON_MODULES_PATH .) #installation root    
	MESSAGE(STATUS "Python modules install path: " "${PYTHON_MODULES_PATH}")
        
ENDMACRO(PIRANHA_PYTHON_SETUP)
