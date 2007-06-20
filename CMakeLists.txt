PROJECT(piranha)
INCLUDE(${CMAKE_ROOT}/Modules/CMakeDependentOption.cmake)
ENABLE_LANGUAGE(C CXX)

# General options
# ---------------------------------------------------------------------------------
OPTION(BUILD_PYRANHA "Build Python bindings." OFF)
CMAKE_DEPENDENT_OPTION(BUILD_PYRANHA_GUI "Build Qt4 gui for pyranha." OFF "BUILD_PYRANHA" OFF)
OPTION(BUILD_MAIN "Build 'main.cpp'." OFF)
OPTION(BUILD_NP "Build numerical series with array lists." ON)
OPTION(BUILD_LNP "Build numerical series with trigonometric lists." OFF)
CMAKE_DEPENDENT_OPTION(BUILD_GSP "Build symbolic series with polynomial coefficients." OFF "BUILD_GMP_SUPPORT" OFF)
CMAKE_DEPENDENT_OPTION(BUILD_TASS "Build TASS17 functions." ON "BUILD_LNP" OFF)

# TODO: remove debugging messages.
MACRO(FIND_BOOST_LIBRARY arg1 arg2)
  IF(WIN32)
    SET(LIBRARY_BASENAME "libboost_${arg2}")
    SET(LIBRARY_NAME_SUFFICES "mgw-mt-d" "mt-d")
  ELSE(WIN32)
    SET(LIBRARY_BASENAME "boost_${arg2}")
    SET(LIBRARY_NAME_SUFFICES "gcc-mt" "mt")
  ENDIF(WIN32)
  SET(BOOST_LIBRARY_NAMES "${LIBRARY_BASENAME}")
  # First we look for plain old name for library.
  FIND_LIBRARY(${arg1} "${LIBRARY_BASENAME}")
  MESSAGE(STATUS "Trying " "${LIBRARY_BASENAME}")
  # Then we iterate over the possible names.
  FOREACH(tmp ${LIBRARY_NAME_SUFFICES})
    SET(BOOST_LIBRARY_NAME "${LIBRARY_BASENAME}-${tmp}")
    # Try to find library only if we did not find it before.
    IF("${arg1}" MATCHES ".*-NOT_FOUND")
      FIND_LIBRARY(${arg1} NAMES "${BOOST_LIBRARY_NAME}")
      MESSAGE(STATUS "Trying " "${BOOST_LIBRARY_NAME}")
    ENDIF("${arg1}" MATCHES ".*-NOT_FOUND")
  ENDFOREACH(tmp ${LIBRARY_NAME_SUFFICES})
  MESSAGE(STATUS "Done macro")
ENDMACRO(FIND_BOOST_LIBRARY)

# Platform switches
IF(WIN32)
  ADD_DEFINITIONS(-D_PIRANHA_WIN32)
  # In MinGW we need the following flags.
  IF(${CMAKE_CXX_COMPILER} MATCHES "(c\\+\\+|g\\+\\+?)")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fdata-sections")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fdata-sections")
    SET(LINK_FLAGS "${LINK_FLAGS} -Wl,--enable-runtime-pseudo-reloc")
  ENDIF(${CMAKE_CXX_COMPILER} MATCHES "(c\\+\\+|g\\+\\+?)")
  SET(LIB_INSTALL_PATH "c:/windows/system32")
ELSE(WIN32)
  SET(LINK_FLAGS "")
  SET(LIB_INSTALL_PATH "lib")
ENDIF(WIN32)

# Compiler flags
# ---------------------------------------------------------------------------------
# Enable verbose reporting of errors when using gcc
SET(GCC_EXTRA_FLAGS "-Wall")
IF(${CMAKE_CXX_COMPILER} MATCHES "c\\+\\+")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_EXTRA_FLAGS}")
ENDIF(${CMAKE_CXX_COMPILER} MATCHES "c\\+\\+")
IF(${CMAKE_C_COMPILER} MATCHES "gcc")
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${GCC_EXTRA_FLAGS}")
ENDIF(${CMAKE_C_COMPILER} MATCHES "gcc")

# Enable further optimization when using icc
# Disable it temporarily, it consumes tons of RAM
#SET(ICC_EXTRA_FLAGS "-ipo")
SET(ICC_EXTRA_FLAGS "")
IF(${CMAKE_CXX_COMPILER} MATCHES "icpc")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ICC_EXTRA_FLAGS}")
ENDIF(${CMAKE_CXX_COMPILER} MATCHES "icpc")
IF(${CMAKE_C_COMPILER} MATCHES "icc")
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${ICC_EXTRA_FLAGS}")
ENDIF(${CMAKE_C_COMPILER} MATCHES "icc")

# Find stuff
# ---------------------------------------------------------------------------------
# Boost general
INCLUDE(${CMAKE_ROOT}/Modules/FindBoost.cmake)
FIND_PACKAGE(Boost REQUIRED)
MESSAGE(STATUS "Boost include path: " "${Boost_INCLUDE_DIRS}")
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
MESSAGE(STATUS "Boost link path: " "${Boost_LIBRARY_DIRS}")
LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})

# Boost filesystem.
FIND_BOOST_LIBRARY(BOOST_FS_LIBRARY "filesystem")
# Boost thread
FIND_BOOST_LIBRARY(BOOST_THREAD_LIBRARY "thread")

# Mandatory linking
SET(MANDATORY_LIBRARIES ${LINK_FLAGS} piranha ${BOOST_FS_LIBRARY} ${BOOST_THREAD_LIBRARY})

# Find GMP and setup support
# Search in this position because it potentially influences many parts of the code -
# i.e., it is not tied to a specific feature.
FIND_LIBRARY(GMP_LIBRARY NAMES "gmp")
FIND_LIBRARY(GMPXX_LIBRARY NAMES "gmpxx")
IF(GMP_LIBRARY AND GMPXX_LIBRARY)
  MESSAGE(STATUS "Found GMP libraries:")
  MESSAGE(STATUS ${GMP_LIBRARY})
  MESSAGE(STATUS ${GMPXX_LIBRARY})
  OPTION(BUILD_GMP_SUPPORT "Build GMP support" ON)
ELSE(GMP_LIBRARY AND GMPXX_LIBRARY)
  MESSAGE(STATUS "GMP libraries not found")
  OPTION(BUILD_GMP_SUPPORT "Build GMP support" OFF)
ENDIF(GMP_LIBRARY AND GMPXX_LIBRARY)
IF(BUILD_GMP_SUPPORT)
  IF(GMP_LIBRARY-NOTFOUND OR GMPXX_LIBRARY-NOTFOUND)
    MESSAGE(FATAL_ERROR "GMP support was requested but libraries were not found.")
  ENDIF(GMP_LIBRARY-NOTFOUND OR GMPXX_LIBRARY-NOTFOUND)
  SET(MANDATORY_LIBRARIES ${MANDATORY_LIBRARIES} ${GMP_LIBRARY} ${GMPXX_LIBRARY})
ENDIF(BUILD_GMP_SUPPORT)

# Fortran
OPTION(BUILD_FORTRAN "Enable support for Fortran routines \(at the moment used by TASS\)" OFF)
IF(BUILD_FORTRAN)
  ENABLE_LANGUAGE(Fortran)
  IF(NOT CMAKE_Fortran_COMPILER)
    MESSAGE(FATAL_ERROR "Fortran support was requested but compiler was not found.")
  ENDIF(NOT CMAKE_Fortran_COMPILER)
  # We need to determine if it is g77 or gfortran: the libraries we will have to link to
  # (manually, since it is fortran inside c++ toolchain) are different. Support for other
  # compilers is a TODO.
  STRING(REGEX MATCH "(gfortran|g77|f77?)" FORTRAN_COMPILER_TYPE ${CMAKE_Fortran_COMPILER})
  IF(${FORTRAN_COMPILER_TYPE} MATCHES "gfortran")
    SET(MANDATORY_LIBRARIES ${MANDATORY_LIBRARIES} gfortran)
  ELSE(${FORTRAN_COMPILER_TYPE} MATCHES "gfortran")
    IF(${FORTRAN_COMPILER_TYPE} MATCHES "g77")
      SET(MANDATORY_LIBRARIES ${MANDATORY_LIBRARIES} g2c)
    ELSE(${FORTRAN_COMPILER_TYPE} MATCHES "g77")
      IF(${FORTRAN_COMPILER_TYPE} MATCHES "f77")
        SET(MANDATORY_LIBRARIES ${MANDATORY_LIBRARIES} g2c)
      ELSE(${FORTRAN_COMPILER_TYPE} MATCHES "f77")
        MESSAGE(FATAL_ERROR "No supported Fortran compiler was found.")
      ENDIF(${FORTRAN_COMPILER_TYPE} MATCHES "f77")
    ENDIF(${FORTRAN_COMPILER_TYPE} MATCHES "g77")
  ENDIF(${FORTRAN_COMPILER_TYPE} MATCHES "gfortran")
  SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -Wall -ffixed-line-length-132")
  ADD_DEFINITIONS(-D_PIRANHA_FORTRAN)
ENDIF(BUILD_FORTRAN)

# TASS
IF(BUILD_TASS)
  ADD_DEFINITIONS(-D_PIRANHA_TASS)
ENDIF(BUILD_TASS)

# Common stuff
# ---------------------------------------------------------------------------------
SET(CMAKE_BUILD_TYPE "Release")
# Find out motion theories files to be installed
SET(THEORIES_INSTALL_PATH "${CMAKE_INSTALL_PREFIX}/share/piranha/theories")
SET(PIRANHA_DEFINITIONS -D_ENABLE_ASSERTS -D_REENTRANT -D_GNU_SOURCE -D_PIRANHA_THEORIES_DIR=\\"${THEORIES_INSTALL_PATH}\\")
ADD_DEFINITIONS(${PIRANHA_DEFINITIONS})

# Add library directory
ADD_SUBDIRECTORY(src)

# Build a main.cpp file included in the source dir
IF(BUILD_MAIN)
  ADD_EXECUTABLE(main main.cpp)
  TARGET_LINK_LIBRARIES(main ${MANDATORY_LIBRARIES})
ENDIF(BUILD_MAIN)

# pyranha: piranha's python bindings
IF(BUILD_PYRANHA)
  ADD_SUBDIRECTORY(pyranha)
ENDIF(BUILD_PYRANHA)

# Install phase
# ---------------------------------------------------------------------------------
# Find and install theories of motion
FILE(GLOB_RECURSE THEORIES theories_of_motion/*.csv theories_of_motion/*.phl theories_of_motion/*.dat)
INSTALL(FILES ${THEORIES} DESTINATION ${THEORIES_INSTALL_PATH})
