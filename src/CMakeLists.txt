# Order is important here: in tass.cpp, for instance, we have destruction of series that depends on
# psymbol_manager to be present. Hence psymbol_manager must go after...
SET(SRC_LIST "")
IF(NOT GCC_CXX_COMPILER AND NOT ICC_CXX_COMPILER AND BUILD_LNP)
  MESSAGE(FATAL_ERROR "The C++ STL extension 'ext/slist' is needed to compile lnp series.")
ELSE(NOT GCC_CXX_COMPILER AND NOT ICC_CXX_COMPILER AND BUILD_LNP)
  MESSAGE(STATUS "Header 'ext/slist' is available.")
ENDIF(NOT GCC_CXX_COMPILER AND NOT ICC_CXX_COMPILER AND BUILD_LNP)
IF(BUILD_TASS)
  IF(BUILD_FORTRAN)
    CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/tass17/posired7.f.template
      ${CMAKE_CURRENT_SOURCE_DIR}/tass17/posired7.f)
    SET(SRC_LIST ${SRC_LIST} tass17/posired7.f)
  ENDIF(BUILD_FORTRAN)
  SET(SRC_LIST ${SRC_LIST} tass17/tass17.cpp)
ENDIF(BUILD_TASS)
IF(BUILD_KEPLERIAN_PROCESSOR)
  SET(SRC_LIST ${SRC_LIST} keplerian_processor/keplerian_processor.cpp)
ENDIF(BUILD_KEPLERIAN_PROCESSOR)
IF(BUILD_SP)
  SET(SRC_LIST ${SRC_LIST} symbol_limiter.cpp)
ENDIF(BUILD_SP)
SET(SRC_LIST ${SRC_LIST} settings_manager.cpp stats.cpp astro.cpp arg_manager.cpp
  stream_manager.cpp psymbol.cpp)
ADD_LIBRARY(piranha SHARED ${SRC_LIST})
IF(BUILD_MAIN)
  ADD_LIBRARY(piranha_static STATIC ${SRC_LIST})
  SET_TARGET_PROPERTIES(piranha_static PROPERTIES LINKER_LANGUAGE CXX)
ENDIF(BUILD_MAIN)
SET_TARGET_PROPERTIES(piranha PROPERTIES LINKER_LANGUAGE CXX)
TARGET_LINK_LIBRARIES(piranha ${MANDATORY_LIBRARIES})
INSTALL(TARGETS piranha
  RUNTIME DESTINATION ${LIB_INSTALL_PATH}
  LIBRARY DESTINATION ${LIB_INSTALL_PATH}
)
