# Copyright (c) 2007, Francesco Biscani, <bluescarni@gmail.com>

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

# Try to find the GSL librairies:
# GSL_FOUND - System has GSL lib
# GSL_INCLUDE_DIR - The GSL include directory
# GSL_LIBRARIES - Libraries needed to use GSL

if (GSL_INCLUDE_DIR AND GSL_LIBRARIES AND GSLCBLAS_LIBRARIES)
  # Already in cache, be silent
  set(GSL_FIND_QUIETLY TRUE)
endif (GSL_INCLUDE_DIR AND GSL_LIBRARIES AND GSLCBLAS_LIBRARIES)

find_path(GSL_INCLUDE_DIR NAMES  gsl)
find_library(GSL_LIBRARIES NAMES gsl)
find_library(GSLCBLAS_LIBRARIES NAMES gslcblas)

if(GSL_INCLUDE_DIR AND GSL_LIBRARIES AND GSLCBLAS_LIBRARIES)
   set(GSL_FOUND TRUE)
endif(GSL_INCLUDE_DIR AND GSL_LIBRARIES AND GSLCBLAS_LIBRARIES)

if(GSL_FOUND)
  if(NOT GSL_FIND_QUIETLY)
    message(STATUS "Found GSL: ${GSL_LIBRARIES}, ${GSLCBLAS_LIBRARIES}")
  endif(NOT GSL_FIND_QUIETLY)
endif(GSL_FOUND)

mark_as_advanced(GSL_INCLUDE_DIR GSL_LIBRARIES GSLCBLAS_LIBRARIES)
