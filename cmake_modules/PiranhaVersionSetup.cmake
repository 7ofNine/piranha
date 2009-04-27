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

# CMake module to setup piranha's version number.

IF(UNIX)
	SET(DATE_BIN "date")
	SET(AWK_BIN "awk")
ELSE(UNIX)
	# Use the utilities from msys.
	SET(UTILS_PATH "c:\\msys\\1.0\\bin")
	SET(DATE_BIN "${UTILS_PATH}\\date.exe")
	SET(AWK_BIN "${UTILS_PATH}\\gawk.exe")
ENDIF(UNIX)

EXECUTE_PROCESS(COMMAND ${DATE_BIN} +%Y RESULT_VARIABLE RETVAL
	OUTPUT_VARIABLE VERSION_MAJOR OUTPUT_STRIP_TRAILING_WHITESPACE)
EXECUTE_PROCESS(COMMAND ${DATE_BIN} +%m RESULT_VARIABLE RETVAL
	OUTPUT_VARIABLE VERSION_MINOR OUTPUT_STRIP_TRAILING_WHITESPACE)
EXECUTE_PROCESS(COMMAND ${DATE_BIN} +%d RESULT_VARIABLE RETVAL
	OUTPUT_VARIABLE VERSION_PATCH OUTPUT_STRIP_TRAILING_WHITESPACE)

MESSAGE(STATUS "Major version number: ${VERSION_MAJOR}")
MESSAGE(STATUS "Minor version number: ${VERSION_MINOR}")
MESSAGE(STATUS "Patch version number: ${VERSION_PATCH}")

SET(PIRANHA_VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")
MESSAGE(STATUS "Piranha version: ${PIRANHA_VERSION}")
