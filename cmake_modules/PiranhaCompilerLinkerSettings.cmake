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

INCLUDE(CheckCXXCompilerFlag)

SET(PIRANHA_EXTRA_LINK_FLAGS "")

IF(MSVC)
        ADD_DEFINITIONS(-DBOOST_ALL_NO_LIB) # don't include boost library versions as #pragma
        ADD_DEFINITIONS(-DNOMINMAX) # get rid of the min/mas macro MS defines
        ADD_DEFINITIONS(-D_SCL_SECURE_NO_WARNINGS)
        ADD_DEFINITIONS(-DWIN32_LEAN_AND_MEAN)  # minimal windows
		ADD_COMPILE_OPTIONS( /std:c++14)
ENDIF()



IF(PIRANHA_EXTRA_LINK_FLAGS)
        MESSAGE(STATUS "Extra link flags: ${PIRANHA_EXTRA_LINK_FLAGS}")
ENDIF(PIRANHA_EXTRA_LINK_FLAGS)

SET(CMAKE_EXE_LINKER_FLAGS    "${CMAKE_EXE_LINKER_FLAGS}    ${PIRANHA_EXTRA_LINK_FLAGS}")
SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${PIRANHA_EXTRA_LINK_FLAGS}")
SET(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} ${PIRANHA_EXTRA_LINK_FLAGS}")
