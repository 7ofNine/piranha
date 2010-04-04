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
IF(CMAKE_COMPILER_IS_GNUCXX)
	MESSAGE(STATUS "Detected GNU compiler")
	# Visibility checks.
	CHECK_CXX_COMPILER_FLAG(-fvisibility-inlines-hidden HAVE_VISIBILITY_INLINES_HIDDEN_FLAG)
	CHECK_CXX_COMPILER_FLAG(-fvisibility=hidden HAVE_VISIBILITY_HIDDEN_FLAG)
	IF(HAVE_VISIBILITY_INLINES_HIDDEN_FLAG AND HAVE_VISIBILITY_HIDDEN_FLAG AND NOT MINGW)
		MESSAGE(STATUS "GCC supports the visibility attributes")
		# SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fvisibility=hidden")
		TRY_COMPILE(WORKING_VISIBILITY_INLINES_HIDDEN ${CMAKE_BINARY_DIR}/compile_tests/ ${CMAKE_SOURCE_DIR}/cmake_modules/simple_main.cpp COMPILE_DEFINITIONS -fvisibility-inlines-hidden)
		IF(WORKING_VISIBILITY_INLINES_HIDDEN)
			MESSAGE(STATUS "GCC version correctly supports the '-fvisibility-inlines-hidden' flag")
			# SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fvisibility-inlines-hidden")
		ELSE(WORKING_VISIBILITY_INLINES_HIDDEN)
			MESSAGE(STATUS "Buggy GCC version, '-fvisibility-inlines-hidden' flag does not work")
		ENDIF(WORKING_VISIBILITY_INLINES_HIDDEN)
	ELSE(HAVE_VISIBILITY_INLINES_HIDDEN_FLAG AND HAVE_VISIBILITY_HIDDEN_FLAG AND NOT MINGW)
		MESSAGE(STATUS "GCC does not support the visibility attributes")
	ENDIF(HAVE_VISIBILITY_INLINES_HIDDEN_FLAG AND HAVE_VISIBILITY_HIDDEN_FLAG AND NOT MINGW)
	# Support for c++0x.
	CHECK_CXX_COMPILER_FLAG(-std=c++0x HAVE_CPP0X_SUPPORT)
	IF(HAVE_CPP0X_SUPPORT)
		MESSAGE(STATUS "Compiler supports c++0x features")
		#ADD_DEFINITIONS(-D_PIRANHA_CPP0X)
		#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
	ENDIF(HAVE_CPP0X_SUPPORT)
	# Extra warnings for the GCC compiler.
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -fmessage-length=0 -Wdisabled-optimization")
	ADD_DEFINITIONS(-D_GNU_SOURCE -D_REENTRANT)
	# Atomic builtins for GCC.
	TRY_COMPILE(WORKING_ATOMIC_BUILTINS ${CMAKE_BINARY_DIR}/compile_tests/ ${CMAKE_SOURCE_DIR}/cmake_modules/gcc_atomic_builtins_test.cpp COMPILE_DEFINITIONS ${CMAKE_CCXX_FLAGS})
	IF(WORKING_ATOMIC_BUILTINS)
		MESSAGE(STATUS "GCC is correctly setup to support atomic builtins.")
		ADD_DEFINITIONS(-D_PIRANHA_GCC_ATOMIC_BUILTINS)
	ELSE(WORKING_ATOMIC_BUILTINS)
		MESSAGE(STATUS "Either this GCC version does not support atomic builtins or the CXXFLAGS are not properly set to enable them.")
		MESSAGE(STATUS "Please note that atomic builtins are supported from version 4.1.0 of GCC and need an appropriate '-march' flag (e.g., at least '-march=i486' on x86 architectures).")
	ENDIF(WORKING_ATOMIC_BUILTINS)
	# MinGW and UNIX specifics.
	IF(MINGW)
		SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fdata-sections -mno-cygwin")
		SET(PIRANHA_EXTRA_LINK_FLAGS "${PIRANHA_EXTRA_LINK_FLAGS} --enable-runtime-pseudo-reloc")
		SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mthreads")
	ELSE(MINGW)
		SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
	ENDIF(MINGW)
ENDIF(CMAKE_COMPILER_IS_GNUCXX)

# Setup link flags.
IF(PIRANHA_EXTRA_LINK_FLAGS)
	MESSAGE(STATUS "Extra link flags: ${PIRANHA_EXTRA_LINK_FLAGS}")
ENDIF(PIRANHA_EXTRA_LINK_FLAGS)
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${PIRANHA_EXTRA_LINK_FLAGS}")
SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${PIRANHA_EXTRA_LINK_FLAGS}")
SET(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} ${PIRANHA_EXTRA_LINK_FLAGS}")
