MACRO(PIRANHA_COMPILER_LINKER_SETTINGS)
  SET(GCC_SPECIFIC_FLAGS "-Wall")
  CHECK_CXX_COMPILER_FLAG(-fvisibility-inlines-hidden __VISIBILITY_INLINES_HIDDEN_FLAG)
  IF(__VISIBILITY_INLINES_HIDDEN_FLAG)
    MESSAGE(STATUS "GCC supports the '-fvisibility-inlines-hidden' flag")
  ELSE(__VISIBILITY_INLINES_HIDDEN_FLAG)
    MESSAGE(STATUS "GCC does not support the '-fvisibility-inlines-hidden' flag")
  ENDIF(__VISIBILITY_INLINES_HIDDEN_FLAG)
  SET(LINK_FLAGS "")
  IF(NOT ${CMAKE_C_COMPILER} MATCHES "gcc" OR NOT ${CMAKE_CXX_COMPILER} MATCHES "(c\\+\\+|g\\+\\+?)")
    MESSAGE(FATAL_ERROR "You need the GCC compiler to build Piranha")
  ENDIF(NOT ${CMAKE_C_COMPILER} MATCHES "gcc" OR NOT ${CMAKE_CXX_COMPILER} MATCHES "(c\\+\\+|g\\+\\+?)")
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${GCC_SPECIFIC_FLAGS}")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_SPECIFIC_FLAGS}")
  IF(ENABLE_HOARD)
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -malign-double")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -malign-double")
  ENDIF(ENABLE_HOARD)
  IF(WIN32)
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fdata-sections")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fdata-sections")
    SET(LINK_FLAGS "${LINK_FLAGS} -Wl,--enable-runtime-pseudo-reloc")
    IF(BUILD_TBB_MULTITHREADING OR ENABLE_HOARD)
      SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mthreads")
      SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mthreads")
    ENDIF(BUILD_TBB_MULTITHREADING OR ENABLE_HOARD)
  ELSE(WIN32)
    IF(BUILD_TBB_MULTITHREADING OR ENABLE_HOARD)
      SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pthread")
      SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
    ENDIF(BUILD_TBB_MULTITHREADING OR ENABLE_HOARD)
  ENDIF(WIN32)
  IF(ENABLE_SSE2)
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -msse2")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse2")
  ENDIF(ENABLE_SSE2)
ENDMACRO(PIRANHA_COMPILER_LINKER_SETTINGS)
