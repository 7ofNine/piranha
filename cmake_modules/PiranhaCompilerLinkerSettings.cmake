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
IF(${CMAKE_CXX_COMPILER} MATCHES "(c\\+\\+|g\\+\\+?)")
	# Visibility checks.
	# TODO: actually check the visibility stuff works, it is not always the case.
	CHECK_CXX_COMPILER_FLAG(-fvisibility-inlines-hidden HAVE_VISIBILITY_INLINES_HIDDEN_FLAG)
	CHECK_CXX_COMPILER_FLAG(-fvisibility=hidden HAVE_VISIBILITY_HIDDEN_FLAG)
	IF(HAVE_VISIBILITY_INLINES_HIDDEN_FLAG AND HAVE_VISIBILITY_HIDDEN_FLAG AND NOT MINGW)
		MESSAGE(STATUS "GCC supports the visibility attributes")
		SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fvisibility=hidden")
		TRY_COMPILE(WORKING_VISIBILITY_INLINES_HIDDEN ${CMAKE_BINARY_DIR}/compile_tests/ ${CMAKE_SOURCE_DIR}/cmake_modules/simple_main.cpp COMPILE_DEFINITIONS -fvisibility-inlines-hidden)
		IF(WORKING_VISIBILITY_INLINES_HIDDEN)
			MESSAGE(STATUS "GCC version correctly supports the '-fvisibility-inlines-hidden' flag")
			SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fvisibility-inlines-hidden")
		ELSE(WORKING_VISIBILITY_INLINES_HIDDEN)
			MESSAGE(STATUS "Buggy GCC version, '-fvisibility-inlines-hidden' flag does not work")
		ENDIF(WORKING_VISIBILITY_INLINES_HIDDEN)
	ELSE(HAVE_VISIBILITY_INLINES_HIDDEN_FLAG AND HAVE_VISIBILITY_HIDDEN_FLAG AND NOT MINGW)
		MESSAGE(STATUS "GCC does not support the visibility attributes")
	ENDIF(HAVE_VISIBILITY_INLINES_HIDDEN_FLAG AND HAVE_VISIBILITY_HIDDEN_FLAG AND NOT MINGW)
	# Support for c++0x.
	CHECK_CXX_COMPILER_FLAG(-std=c++0x HAVE_CPP0X_SUPPORT)
	IF(HAVE_CPP0X_SUPPORT)
		#MESSAGE(STATUS "Enabling support for c++0x features.")
		#ADD_DEFINITIONS(-D_PIRANHA_CPP0X)
		#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
	ENDIF(HAVE_CPP0X_SUPPORT)
	# Extra warnings for the GCC compiler.
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -fmessage-length=0 -Wdisabled-optimization")
	ADD_DEFINITIONS(-D_GNU_SOURCE)
	IF(BUILD_MULTITHREADING)
		SET(PIRANHA_DEFINITIONS ${PIRANHA_DEFINITIONS} -D_REENTRANT)
	ENDIF(BUILD_MULTITHREADING)
	# MinGW and UNIX specifics.
	IF(WIN32)
		SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fdata-sections -mno-cygwin")
		SET(PIRANHA_EXTRA_LINK_FLAGS "${PIRANHA_EXTRA_LINK_FLAGS} --enable-runtime-pseudo-reloc")
		IF(BUILD_MULTITHREADING)
			SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mthreads")
		ENDIF(BUILD_MULTITHREADING)
	ELSE(WIN32)
		IF(BUILD_MULTITHREADING)
			SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
		ENDIF(BUILD_MULTITHREADING)
	ENDIF(WIN32)
	SET(GNU_COMPILER TRUE)
ENDIF(${CMAKE_CXX_COMPILER} MATCHES "(c\\+\\+|g\\+\\+?)")

# Setup link flags.
IF(PIRANHA_EXTRA_LINK_FLAGS)
	MESSAGE(STATUS "Extra link flags: ${PIRANHA_EXTRA_LINK_FLAGS}")
ENDIF(PIRANHA_EXTRA_LINK_FLAGS)
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${PIRANHA_EXTRA_LINK_FLAGS}")
SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${PIRANHA_EXTRA_LINK_FLAGS}")
SET(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} ${PIRANHA_EXTRA_LINK_FLAGS}")
