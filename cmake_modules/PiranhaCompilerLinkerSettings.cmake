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

MACRO(PIRANHA_COMPILER_LINKER_SETTINGS)
  SET(GCC_SPECIFIC_FLAGS "-Wall -Wextra -fmessage-length=0 -Wdisabled-optimization")
  # Visibility checks.
  CHECK_CXX_COMPILER_FLAG(-fvisibility-inlines-hidden __VISIBILITY_INLINES_HIDDEN_FLAG)
  CHECK_CXX_COMPILER_FLAG(-fvisibility=hidden __VISIBILITY_HIDDEN_FLAG)
  IF(__VISIBILITY_INLINES_HIDDEN_FLAG AND __VISIBILITY_HIDDEN_FLAG)
    MESSAGE(STATUS "GCC supports the visibility attributes")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fvisibility-inlines-hidden -fvisibility=hidden")
  ELSE(__VISIBILITY_INLINES_HIDDEN_FLAG AND __VISIBILITY_HIDDEN_FLAG)
    MESSAGE(STATUS "GCC does not support the visibility attributes")
  ENDIF(__VISIBILITY_INLINES_HIDDEN_FLAG AND __VISIBILITY_HIDDEN_FLAG)
  SET(LINK_FLAGS "")
  IF(${CMAKE_CXX_COMPILER} MATCHES "(c\\+\\+|g\\+\\+?)")
    SET(GNU_COMPILER TRUE)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_SPECIFIC_FLAGS}")
  ENDIF(${CMAKE_CXX_COMPILER} MATCHES "(c\\+\\+|g\\+\\+?)")
  IF(ENABLE_HOARD AND NOT ${POINTER_SIZE} EQUAL 8)
    IF(GNU_COMPILER)
      SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -malign-double")
      SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -malign-double")
    ENDIF(GNU_COMPILER)
  ENDIF(ENABLE_HOARD AND NOT ${POINTER_SIZE} EQUAL 8)
  IF(MINGW)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fdata-sections")
    SET(LINK_FLAGS "${LINK_FLAGS} -Wl,--enable-runtime-pseudo-reloc")
    IF(BUILD_TBB_MULTITHREADING OR ENABLE_HOARD)
      SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mthreads")
    ENDIF(BUILD_TBB_MULTITHREADING OR ENABLE_HOARD)
  ELSE(MINGW)
    IF(GNU_COMPILER)
      IF(BUILD_TBB_MULTITHREADING OR ENABLE_HOARD)
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
      ENDIF(BUILD_TBB_MULTITHREADING OR ENABLE_HOARD)
    ENDIF(GNU_COMPILER)
  ENDIF(MINGW)
  IF(ENABLE_SSE2)
    IF(GNU_COMPILER)
      SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse2")
    ENDIF(GNU_COMPILER)
  ENDIF(ENABLE_SSE2)
ENDMACRO(PIRANHA_COMPILER_LINKER_SETTINGS)
