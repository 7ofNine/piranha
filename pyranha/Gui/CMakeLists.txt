
message(STATUS "module path ${CMAKE_MODULE_PATH}")
INCLUDE(FindPySide2)

FIND_PACKAGE(PySide2 REQUIRED)

SET(ui_list main_window.ui)
#SET(rc_list __resources.qrc)
PYSIDE2_WRAP_UI(pyranha_ui ${ui_list})
#PYSIDE2_WRAP_RC(pyranha_rc ${rc_list})

SET(pyranha_py main_window.py series_draw_area.py)

MESSAGE(STATUS "pyranha_ui is " ${pyranha_ui})
#MESSAGE(STATUS "pyranha_rc is " ${pyranha_rc})

#SET(IPYTHON_GUI_TOOLKIT "-q4thread")
SET(IPYTHON_GUI_TOOLKIT "--gui=qt")
CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/windows_package/template.console.xml ${CMAKE_SOURCE_DIR}/windows_package/console.xml)

#INSTALL(FILES __init__.py panel_widget.py ${pyranha_ui} ${pyranha_rc} DESTINATION ${PYRANHA_INSTALL_PATH}/Gui/)
INSTALL(FILES __init__.py ${pyranha_py} ${pyranha_ui} DESTINATION ${PYRANHA_INSTALL_PATH}/Gui/)
