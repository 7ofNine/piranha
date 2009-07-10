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

EXECUTE_PROCESS(COMMAND git log --no-color -n1 --date=short --pretty=format:%ad OUTPUT_VARIABLE PIRANHA_GIT_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE)
EXECUTE_PROCESS(COMMAND git log --no-color -n1 --pretty=format:%h OUTPUT_VARIABLE PIRANHA_GIT_REVISION OUTPUT_STRIP_TRAILING_WHITESPACE)

STRING(REGEX REPLACE "-" "." PIRANHA_VERSION ${PIRANHA_GIT_VERSION})

MESSAGE(STATUS "Piranha version: ${PIRANHA_VERSION}")
