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
