MACRO(PIRANHA_COMPILER_LINKER_SETTINGS)
  # Setup variables.
  SET(ICC_C_COMPILER FALSE)
  SET(GCC_C_COMPILER FALSE)
  SET(ICC_CXX_COMPILER FALSE)
  SET(GCC_CXX_COMPILER FALSE)
  SET(GCC_SPECIFIC_FLAGS "-Wall")
  # Disable it temporarily, it consumes tons of RAM.
  #SET(ICC_SPECIFIC_FLAGS "-ipo")
  SET(ICC_SPECIFIC_FLAGS "")
  SET(LINK_FLAGS "")
  # Establish compiler type. Supported: GCC, ICC.
  IF(${CMAKE_C_COMPILER} MATCHES "gcc")
    SET(GCC_C_COMPILER TRUE)
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${GCC_SPECIFIC_FLAGS}")
    IF(WIN32)
      SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fdata-sections")
      IF(BUILD_TBB_MULTITHREADING)
        SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mthreads")
      ENDIF(BUILD_TBB_MULTITHREADING)
    ELSE(WIN32)
      IF(BUILD_TBB_MULTITHREADING)
        SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pthread")
      ENDIF(BUILD_TBB_MULTITHREADING)
    ENDIF(WIN32)
  ENDIF(${CMAKE_C_COMPILER} MATCHES "gcc")
  IF(${CMAKE_CXX_COMPILER} MATCHES "(c\\+\\+|g\\+\\+?)")
    SET(GCC_CXX_COMPILER TRUE)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_SPECIFIC_FLAGS}")
    IF(WIN32)
      SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fdata-sections")
      SET(LINK_FLAGS "${LINK_FLAGS} -Wl,--enable-runtime-pseudo-reloc")
      IF(BUILD_TBB_MULTITHREADING)
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mthreads")
      ENDIF(BUILD_TBB_MULTITHREADING)
    ELSE(WIN32)
      IF(BUILD_TBB_MULTITHREADING)
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
      ENDIF(BUILD_TBB_MULTITHREADING)
    ENDIF(WIN32)
  ENDIF(${CMAKE_CXX_COMPILER} MATCHES "(c\\+\\+|g\\+\\+?)")
  IF(${CMAKE_C_COMPILER} MATCHES "icc")
    SET(ICC_C_COMPILER TRUE)
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${ICC_SPECIFIC_FLAGS}")
  ENDIF(${CMAKE_C_COMPILER} MATCHES "icc")
  IF(${CMAKE_CXX_COMPILER} MATCHES "icpc")
    SET(ICC_CXX_COMPILER TRUE)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ICC_SPECIFIC_FLAGS}")
  ENDIF(${CMAKE_CXX_COMPILER} MATCHES "icpc")
ENDMACRO(PIRANHA_COMPILER_LINKER_SETTINGS)
