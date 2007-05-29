find_package(Qt4 REQUIRED)

FIND_PROGRAM(PYUIC4BINARY pyuic4)
IF(${PYUIC4BINARY} MATCHES "PYUIC4BINARY-NOTFOUND")
  MESSAGE(FATAL_ERROR "pyuic4 not found, please install PyQt4 or disable gui.")
ELSE(${PYUIC4BINARY} MATCHES "PYUIC4BINARY-NOTFOUND")
  MESSAGE(STATUS "pyuic4 binary is located at " ${PYUIC4BINARY})
ENDIF(${PYUIC4BINARY} MATCHES "PYUIC4BINARY-NOTFOUND")

FIND_PROGRAM(PYRCC4BINARY pyrcc4)
IF(${PYRCC4BINARY} MATCHES "PYRCC4BINARY-NOTFOUND")
  MESSAGE(FATAL_ERROR "pyrcc4 not found, please install PyQt4 or disable gui.")
ELSE(${PYRCC4BINARY} MATCHES "PYRCC4BINARY-NOTFOUND")
  MESSAGE(STATUS "pyrcc4 binary is located at " ${PYRCC4BINARY})
ENDIF(${PYRCC4BINARY} MATCHES "PYRCC4BINARY-NOTFOUND")

MACRO(PYQT4_WRAP_UI outfiles)
  FOREACH(it ${ARGN})
    GET_FILENAME_COMPONENT(outfile ${it} NAME_WE)
    GET_FILENAME_COMPONENT(infile ${it} ABSOLUTE)
    SET(outfile ${CMAKE_CURRENT_BINARY_DIR}/ui_${outfile}.py)
    ADD_CUSTOM_TARGET(${it} ALL
      DEPENDS ${outfile}
    )
    ADD_CUSTOM_COMMAND(OUTPUT ${outfile}
      COMMAND ${PYUIC4BINARY} ${infile} -o ${outfile}
      MAIN_DEPENDENCY ${infile}
    )
    SET(${outfiles} ${${outfiles}} ${outfile})
  ENDFOREACH(it)
ENDMACRO (PYQT4_WRAP_UI)

MACRO(PYQT4_WRAP_RC outfiles)
  FOREACH(it ${ARGN})
    GET_FILENAME_COMPONENT(outfile ${it} NAME_WE)
    GET_FILENAME_COMPONENT(infile ${it} ABSOLUTE)
    SET(outfile ${CMAKE_CURRENT_BINARY_DIR}/${outfile}_rc.py)
    ADD_CUSTOM_TARGET(${it} ALL
      DEPENDS ${outfile}
    )
    ADD_CUSTOM_COMMAND(OUTPUT ${outfile}
      COMMAND ${PYRCC4BINARY} ${infile} -o ${outfile}
      MAIN_DEPENDENCY ${infile}
    )
    SET(${outfiles} ${${outfiles}} ${outfile})
  ENDFOREACH(it)
ENDMACRO (PYQT4_WRAP_RC)

SET(ui_list __panelWidget.ui)
SET(rc_list __resources.qrc)
PYQT4_WRAP_UI(pyranha_ui ${ui_list})
PYQT4_WRAP_RC(pyranha_rc ${rc_list})

MESSAGE(STATUS "pyranha_ui is " ${pyranha_ui})
MESSAGE(STATUS "pyranha_rc is " ${pyranha_rc})

INSTALL(FILES __init__.py __panelWidget.py ${pyranha_ui} ${pyranha_rc} DESTINATION ${PYRANHA_INSTALL_PATH}/Gui/)
