# Copyright (C) 2007, 2008 by Francesco Biscani; 2017 H=by Hartmuth Gutsche
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

# CMake module to setup piranha's version number.

FIND_PROGRAM(PIRANHA_GIT_EXECUTABLE git PATHS "C:\\Program Files (x86)\\" "C:\\Program Files\\" PATH_SUFFIXES "git" "git\\bin" DOC "Path to the git binary.")

IF(PIRANHA_GIT_EXECUTABLE)
        MESSAGE(STATUS "Git executable: ${PIRANHA_GIT_EXECUTABLE}")

        # Version number setup.
        SET(PIRANHA_GIT_ARGS "log" "--no-color" "-n1" "--date=short" "--pretty=format:%ad")
        EXECUTE_PROCESS(COMMAND ${PIRANHA_GIT_EXECUTABLE} ${PIRANHA_GIT_ARGS} WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} OUTPUT_VARIABLE PIRANHA_GIT_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE)
#        STRING(REGEX REPLACE "-" "." PIRANHA_VERSION ${PIRANHA_GIT_VERSION})
        MESSAGE(STATUS "Piranha version: ${PIRANHA_GIT_VERSION}")
        SET(PIRANHA_VERSION_DEFINE "#define PIRANHA_VERSION ${PIRANHA_GIT_VERSION}")
        SET(PIRANHA_VERSION_STRING "${PIRANHA_GIT_VERSION}")
		
		SET(PIRANHA_GIT_ARGS "log" "--no-color" "-n1" "--date=short" "--pretty=format:%H")
		EXECUTE_PROCESS(COMMAND ${PIRANHA_GIT_EXECUTABLE} ${PIRANHA_GIT_ARGS} WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} OUTPUT_VARIABLE PIRANHA_GIT_REVISION OUTPUT_STRIP_TRAILING_WHITESPACE)
		MESSAGE(STATUS "Piranha revision: ${PIRANHA_GIT_REVISION}") 
		
ELSE(PIRANHA_GIT_EXECUTABLE)
        MESSAGE(STATUS "Git executable: not found. Can't determine Piranha version")
        SET(PIRANHA_VERSION_STRING "Undetermined")
		SET(PIRANHA_GIT_REVISION 0)
ENDIF(PIRANHA_GIT_EXECUTABLE)
