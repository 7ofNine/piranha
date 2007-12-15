# Fortran
MACRO(PIRANHA_FORTRAN_SETUP)
  ENABLE_LANGUAGE(Fortran)
  IF(NOT CMAKE_Fortran_COMPILER)
    MESSAGE(FATAL_ERROR "Fortran support was requested but compiler was not found.")
  ENDIF(NOT CMAKE_Fortran_COMPILER)
  # We need to determine the compiler: the libraries we will have to link to
  # (manually, since it is fortran inside c++ toolchain) are different. Support for other
  # compilers can be added.
  STRING(REGEX MATCH "(gfortran|g77|f77?)" FORTRAN_COMPILER_TYPE ${CMAKE_Fortran_COMPILER})
  IF(${FORTRAN_COMPILER_TYPE} MATCHES "gfortran")
    # This seems a bit of a hack, but we cannot use find_library because libgfortran, if existing,
    # is not in standard paths.
    CHECK_CXX_COMPILER_FLAG(-lgfortran __LIBGFORTRAN_FLAG)
    IF(__LIBGFORTRAN_FLAG)
      MESSAGE(STATUS "libgfortran was detected, will link against it.")
      SET(MANDATORY_LIBRARIES ${MANDATORY_LIBRARIES} gfortran)
    ENDIF(__LIBGFORTRAN_FLAG)
    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O2 -Wall -ffixed-line-length-132")
  ELSE(${FORTRAN_COMPILER_TYPE} MATCHES "gfortran")
    IF(${FORTRAN_COMPILER_TYPE} MATCHES "g77")
      SET(MANDATORY_LIBRARIES ${MANDATORY_LIBRARIES} g2c)
      SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O2 -Wall -ffixed-line-length-132")
    ELSE(${FORTRAN_COMPILER_TYPE} MATCHES "g77")
      IF(${FORTRAN_COMPILER_TYPE} MATCHES "f77")
        SET(MANDATORY_LIBRARIES ${MANDATORY_LIBRARIES} g2c)
        SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O2 -Wall -ffixed-line-length-132")
      ELSE(${FORTRAN_COMPILER_TYPE} MATCHES "f77")
        MESSAGE(FATAL_ERROR "No supported Fortran compiler was found.")
      ENDIF(${FORTRAN_COMPILER_TYPE} MATCHES "f77")
    ENDIF(${FORTRAN_COMPILER_TYPE} MATCHES "g77")
  ENDIF(${FORTRAN_COMPILER_TYPE} MATCHES "gfortran")
  ADD_DEFINITIONS(-D_PIRANHA_FORTRAN)
ENDMACRO(PIRANHA_FORTRAN_SETUP)
