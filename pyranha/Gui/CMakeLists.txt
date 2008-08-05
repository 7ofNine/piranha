INCLUDE(FindPyQt4)

FIND_PACKAGE(PyQt4 REQUIRED)

SET(ui_list main_window.ui)
#SET(rc_list __resources.qrc)
PYQT4_WRAP_UI(pyranha_ui ${ui_list})
#PYQT4_WRAP_RC(pyranha_rc ${rc_list})

SET(pyranha_py main_window.py series_draw_area.py)

MESSAGE(STATUS "pyranha_ui is " ${pyranha_ui})
MESSAGE(STATUS "pyranha_rc is " ${pyranha_rc})

SET(IPYTHON_GUI_TOOLKIT "-q4thread")
CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/package/template.console.xml ${CMAKE_SOURCE_DIR}/package/console.xml)

#INSTALL(FILES __init__.py panel_widget.py ${pyranha_ui} ${pyranha_rc} DESTINATION ${PYRANHA_INSTALL_PATH}/Gui/)
INSTALL(FILES __init__.py ${pyranha_py} ${pyranha_ui} DESTINATION ${PYRANHA_INSTALL_PATH}/Gui/)
