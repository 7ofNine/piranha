# Order is important here: in tass.cpp, for instance, we have destruction of series that depends on
# psymbol_manager to be present. Hence psymbol_manager must go after...
SET(SRC_LIST "")
IF(BUILD_TASS)
  IF(BUILD_FORTRAN)
    CONFIGURE_FILE(./tass17/posired7.f.template src/tass17/posired7.f)
    SET(SRC_LIST ${SRC_LIST} tass17/posired7.f)
  ENDIF(BUILD_FORTRAN)
  SET(SRC_LIST ${SRC_LIST} tass17/tass17.cpp)
ENDIF(BUILD_TASS)
SET(SRC_LIST ${SRC_LIST} settings_manager.cpp stats.cpp psymbol.cpp astro.cpp arg_manager.cpp
  stream_manager.cpp psymbol_manager.cpp)
ADD_LIBRARY(piranha SHARED ${SRC_LIST})
SET_TARGET_PROPERTIES(piranha PROPERTIES LINKER_LANGUAGE CXX)
TARGET_LINK_LIBRARIES(piranha ${MANDATORY_LIBRARIES})
INSTALL(TARGETS piranha
  RUNTIME DESTINATION ${LIB_INSTALL_PATH}
  LIBRARY DESTINATION ${LIB_INSTALL_PATH}
)
