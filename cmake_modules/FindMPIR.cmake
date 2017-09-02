# Originally copied from the KDE project repository:
# http://websvn.kde.org/trunk/KDE/kdeutils/cmake/modules/FindGMP.cmake?view=markup&pathrev=675218

# Copyright (c) 2006, Laurent Montel, <montel@kde.org>
# Copyright (c) 2007, 2008 Francesco Biscani, <bluescarni@gmail.com>
# COpyright (c) 2017 Hartmuth Gutsche

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
# 1. Redistributions of source code must retain the copyright
#    notice, this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
# 3. The name of the author may not be used to endorse or promote products 
#    derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
# IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
# OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
# IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
# NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
# THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# ------------------------------------------------------------------------------------------


# Try to find the MPIR librairies:
# MPIR_FOUND        - System has MPIR lib
# MPIR_INCLUDE_DIR  - The MPIR include directory
# MPIR_LIBRARIES    - Libraries needed to use MPIR and its C++ interface
#
# for MSVC we use the MPIR, and windows MPFR implementations of GMP/MPFR
# see git@github.com:wbhart/mpir.git and https://github.com/BrianGladman/mpfr.git
# i.e. all names are changed to MPIR, originally this was Gnu GMP
# 


IF(MPIR_INCLUDE_DIR AND MPIR_LIBRARIES)
	# Already in cache, be silent
	SET(MPIR_FIND_QUIETLY TRUE)
ENDIF(MPIR_INCLUDE_DIR AND MPIR_LIBRARIES)

# We are only looking for the DLL version. The dll version should also include the C++ interface
# It gets linked via the *.lib stub
# we also only want the 64 Bit versions
FIND_PATH(MPIR_INCLUDE_DIR NAMES mpir.h PATHS D:\\Dev\\mpir\\dll\\x64\\Release) # we don't use this file but thats where gmp.h, gmpxx.h are, too
FIND_LIBRARY(MPIR_LIBRARIES NAMES mpir.lib PATHS D:\\Dev\\mpir\\dll\\x64\\Release)

IF(MPIR_INCLUDE_DIR AND MPIR_LIBRARIES)
	SET(MPIR_FOUND TRUE)
ENDIF(MPIR_INCLUDE_DIR AND MPIR_LIBRARIES)

IF(MPIR_FOUND)
	IF(NOT MPIR_FIND_QUIETLY)
		MESSAGE(STATUS "Found MPIR: ${MPIR_LIBRARIES}")
	ENDIF(NOT MPIR_FIND_QUIETLY)
ELSE(MPIR_FOUND)
	IF(MPIR_FIND_REQUIRED)
		MESSAGE(SEND_ERROR "Could NOT find MPIR")
	ENDIF(MPIR_FIND_REQUIRED)
ENDIF(MPIR_FOUND)

MARK_AS_ADVANCED(MPIR_INCLUDE_DIR MPIR_LIBRARIES)
