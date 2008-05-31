INCLUDE(FindPyQt4)

find_package(Qt4 REQUIRED)
find_package(PyQt4 REQUIRED)

SET(ui_list panel_widget.ui)
#SET(rc_list __resources.qrc)
PYQT4_WRAP_UI(pyranha_ui ${ui_list})
#PYQT4_WRAP_RC(pyranha_rc ${rc_list})

MESSAGE(STATUS "pyranha_ui is " ${pyranha_ui})
MESSAGE(STATUS "pyranha_rc is " ${pyranha_rc})

#INSTALL(FILES __init__.py panel_widget.py ${pyranha_ui} ${pyranha_rc} DESTINATION ${PYRANHA_INSTALL_PATH}/Gui/)
INSTALL(FILES __init__.py panel_widget.py ${pyranha_ui} DESTINATION ${PYRANHA_INSTALL_PATH}/Gui/)