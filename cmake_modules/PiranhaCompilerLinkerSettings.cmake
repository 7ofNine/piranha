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
INCLUDE(CheckTypeSize)

MACRO(PIRANHA_COMPILER_LINKER_SETTINGS)
	SET(PIRANHA_EXTRA_LINK_FLAGS "")
	CHECK_TYPE_SIZE("void *" POINTER_SIZE)
	MESSAGE(STATUS "Pointer size is " ${POINTER_SIZE})
	IF(${CMAKE_CXX_COMPILER} MATCHES "(c\\+\\+|g\\+\\+?)")
		# Visibility checks.
		CHECK_CXX_COMPILER_FLAG(-fvisibility-inlines-hidden HAVE_VISIBILITY_INLINES_HIDDEN_FLAG)
		CHECK_CXX_COMPILER_FLAG(-fvisibility=hidden HAVE_VISIBILITY_HIDDEN_FLAG)
		IF(HAVE_VISIBILITY_INLINES_HIDDEN_FLAG AND HAVE_VISIBILITY_HIDDEN_FLAG)
			MESSAGE(STATUS "GCC supports the visibility attributes")
			SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fvisibility-inlines-hidden -fvisibility=hidden")
		ELSE(HAVE_VISIBILITY_INLINES_HIDDEN_FLAG AND HAVE_VISIBILITY_HIDDEN_FLAG)
			MESSAGE(STATUS "GCC does not support the visibility attributes")
		ENDIF(HAVE_VISIBILITY_INLINES_HIDDEN_FLAG AND HAVE_VISIBILITY_HIDDEN_FLAG)
		CHECK_CXX_COMPILER_FLAG(-std=c++0x HAVE_CPP0X_SUPPORT)
		IF(HAVE_CPP0X_SUPPORT)
			#MESSAGE(STATUS "Enabling support for c++0x features.")
			#ADD_DEFINITIONS(-D_PIRANHA_CPP0X)
			#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
		ENDIF(HAVE_CPP0X_SUPPORT)
		# Extra warnings for the GCC compiler.
		SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -fmessage-length=0 -Wdisabled-optimization")
		ADD_DEFINITIONS(-D_GNU_SOURCE)
		SET(GNU_COMPILER TRUE)
	ENDIF(${CMAKE_CXX_COMPILER} MATCHES "(c\\+\\+|g\\+\\+?)")
	# Hoard requires double variables to be aligned on two-word boundaries.
	# On 64-bit archs this should be the default.
	IF(ENABLE_HOARD AND NOT ${POINTER_SIZE} EQUAL 8)
		IF(GNU_COMPILER)
			SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -malign-double")
			SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -malign-double")
		ENDIF(GNU_COMPILER)
	ENDIF(ENABLE_HOARD AND NOT ${POINTER_SIZE} EQUAL 8)
	# Some MinGW setup.
	IF(MINGW)
		SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fdata-sections")
		SET(PIRANHA_EXTRA_LINK_FLAGS "${PIRANHA_EXTRA_LINK_FLAGS} --enable-runtime-pseudo-reloc")
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
	# SSE2 support.
	IF(ENABLE_SSE2)
		IF(GNU_COMPILER)
			CHECK_CXX_COMPILER_FLAG(-msse2 HAVE_SSE2_FLAG)
			IF(HAVE_SSE2_FLAG)
				SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse2")
			ELSE(HAVE_SSE2_FLAG)
				MESSAGE(FATAL_ERROR "SSE2 support requested but compiler does not support the '-msse2' flag.")
			ENDIF(HAVE_SSE2_FLAG)
		ENDIF(GNU_COMPILER)
	ENDIF(ENABLE_SSE2)
	# Setup link flags.
	IF(PIRANHA_EXTRA_LINK_FLAGS)
		MESSAGE(STATUS "Extra link flags: ${PIRANHA_EXTRA_LINK_FLAGS}")
	ENDIF(PIRANHA_EXTRA_LINK_FLAGS)
	SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${PIRANHA_EXTRA_LINK_FLAGS}")
	SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${PIRANHA_EXTRA_LINK_FLAGS}")
	SET(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} ${PIRANHA_EXTRA_LINK_FLAGS}")
ENDMACRO(PIRANHA_COMPILER_LINKER_SETTINGS)
