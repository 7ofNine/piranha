# Copyright (C) 2007, 2008 by Francesco Biscani; 2017 Hartmuth Gutsche
#               2021 Hartmuth Gutsche
#
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

MACRO(PIRANHA_BOOST_SETUP)

    # Boost find parameters
    SET(Boost_USE_RELEASE_LIBS ON)
    SET(Boost_USE_STATIC_LIBS OFF)
    SET(Boost_USE_STATIC_RUNTIME OFF)
    set(Boost_USE_MULTITHREADED ON)
    set(Boost_USE_DEBUG_PYTHON ON)
    
    If(BUILD_PYRANHA)
        set(PYRANHA_BOOST_PYTHON_TARGET ${PIRANHA_PYTHON_BASE})
    endif()
    
    FIND_PACKAGE(Boost 1.74.0 REQUIRED COMPONENTS thread ${PYRANHA_BOOST_PYTHON_TARGET})
    
    if(Boost_FOUND)
        MESSAGE(STATUS "Found Boost libraries: ${Boost_LIBRARIES}")
        SET(PYRANHA_BOOST_LIBRARIES "${Boost_LIBRARIES}")  #push them up in scope??
        MESSAGE(STATUS "thread library: ${Boost_thread_LIBRARY_RELEASE}")
    else()
        MESSAGE(FATAL "Couldn't find Boost libraries")
    endif()
ENDMACRO(PIRANHA_BOOST_SETUP)
