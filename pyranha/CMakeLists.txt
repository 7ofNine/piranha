#
# pyranha/CMakeList
#
SET(PYRANHA_INSTALL_PATH ${PYTHON_MODULES_PATH}/pyranha)
MESSAGE (STATUS "pyranha install path ${PYRANHA_INSTALL_PATH}")

ADD_COMPILE_OPTIONS(/bigobj)
ADD_COMPILE_OPTIONS(/wd4244) # Boost python shows a lot of size warnings with 64 bit compilation. Boost should fix that
ADD_COMPILE_OPTIONS(/wd4251)
ADD_COMPILE_DEFINITIONS(PIRANHA_DLL_IMPORT_API BOOST_BIND_GLOBAL_PLACEHOLDERS)  # here we import piranha.dll until boost/python fixes the deprecated warnings
INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIR})  # everything under here needs the python declarations
link_directories(${PIRANHA_PYTHON_LIBRARIES_DIR})    # pyconfig.h sets required pylib (or pylib_d) which gets messed up by boost wrap_python.h when _DEBUG is set, it temporary unsets it. This requires the python library dir explicitly available
                                                     # the pyconfig.h only contains the names not the path!!! We don't set _DEBUG hence we have to provide the library dir separately


MACRO(PYTHON_CPPMODULE arg1)
    SET_TARGET_PROPERTIES(${arg1} PROPERTIES PREFIX "")
    SET_TARGET_PROPERTIES(${arg1} PROPERTIES SUFFIX ".pyd")
    SET_TARGET_PROPERTIES(${arg1} PROPERTIES FOLDER "Pyranha") # we want to put these projects into a Pyranha folder in VS
    TARGET_LINK_LIBRARIES(${arg1} ${MANDATORY_LIBRARIES} ${PIRANHA_PYTHON_LIBRARIES} ${PYRANHA_BOOST_LIBRARIES} ${ARGN})
    TARGET_INCLUDE_DIRECTORIES(${arg1} PRIVATE ${PIRANHA_PYTHON_INCLUDE_DIR})
    STRING(REGEX REPLACE "^_" "" DIR_ARG ${arg1})
#   MESSAGE(STATUS "Install files: ${arg1}")
    INSTALL(TARGETS ${arg1} RUNTIME DESTINATION ${PYRANHA_INSTALL_PATH}/${DIR_ARG} )
ENDMACRO(PYTHON_CPPMODULE)


MACRO(PYRANHA_MANIPULATOR_INIT arg1)
    STRING(REGEX REPLACE "^_" "" DIR_ARG ${arg1})
    SET(MODULE_NAME ${DIR_ARG})  # we use DIR_ARG because it is the name without "_". Needed in the __init__.py file
    CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/../template.__init__.py ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py)
    INSTALL(FILES __init__.py DESTINATION ${PYRANHA_INSTALL_PATH}/${DIR_ARG})
ENDMACRO(PYRANHA_MANIPULATOR_INIT)


MACRO(PYRANHA_MANIPULATOR arg1)
    PYTHON_CPPMODULE(${arg1})
    PYRANHA_MANIPULATOR_INIT(${arg1})
ENDMACRO(PYRANHA_MANIPULATOR)

SET(PYRANHA_HEADERS
    args_tuple.h
    boost_python_container_conversions.h
    boost_python_p_exceptions.h
    commons.h
    exceptions.h
    mp_classes.h
    series_instantiations.h
    stl_containers.h)

SET(PYRANHA_PYTHON_COMMON
    detail.py
    global.__init__.py
    template.__init__.py
    template.global.__init__.py
    )
    
SOURCE_GROUP(TREE ${CMAKE_CURRENT_SOURCE_DIR} PREFIX "Header Files" FILES ${PYRANHA_HEADERS}) 
SOURCE_GROUP(TREE ${CMAKE_CURRENT_SOURCE_DIR} PREFIX "Common Python Files" FILES ${PYRANHA_PYTHON_COMMON}) 
ADD_LIBRARY(Common STATIC ${PYRANHA_HEADERS} ${PYRANHA_PYTHON_COMMON}) 
SET_TARGET_PROPERTIES(Common PROPERTIES LINKER_LANGUAGE CXX) #Hack to make cmake error go away
SET_TARGET_PROPERTIES(Common PROPERTIES FOLDER "Pyranha") # Put Common under Pyranha Folder
 

#MESSAGE(STATUS ${CMAKE_SOURCE_DIR} )
# Mandatory targets.
ADD_SUBDIRECTORY(Core)
SET(MODULE_LIST \ \"Core\",\ )

ADD_SUBDIRECTORY(Math)
SET(MODULE_LIST ${MODULE_LIST}"Math",\ )

ADD_SUBDIRECTORY(Celmec)
SET(MODULE_LIST ${MODULE_LIST}"Celmec",\ )

ADD_SUBDIRECTORY(Test)
SET(MODULE_LIST ${MODULE_LIST}"Test",\ )

ADD_SUBDIRECTORY(Truncators)
SET(MODULE_LIST ${MODULE_LIST}"Truncators",\ )

# GUI (optional). There is not much in it.
IF(BUILD_PYRANHA_GUI)
    SET(MODULE_LIST ${MODULE_LIST}"Gui",\ )
    ADD_SUBDIRECTORY(Gui)
ENDIF(BUILD_PYRANHA_GUI)

# Manipulators.
IF(BUILD_PYRANHA_DFS)
    SET(MODULE_LIST ${MODULE_LIST}"Dfs",\ )
    SET(MANIPULATOR_LIST ${MANIPULATOR_LIST}"Dfs",\ )
    ADD_SUBDIRECTORY(Dfs)
ENDIF(BUILD_PYRANHA_DFS)

IF(BUILD_PYRANHA_DPOLY)
    SET(MODULE_LIST ${MODULE_LIST}"Dpoly",\ )
    SET(MANIPULATOR_LIST ${MANIPULATOR_LIST}"Dpoly",\ )
    ADD_SUBDIRECTORY(Dpoly)
ENDIF(BUILD_PYRANHA_DPOLY)

IF(BUILD_PYRANHA_ZPOLY)
    SET(MODULE_LIST ${MODULE_LIST}"Zpoly",\ )
    SET(MANIPULATOR_LIST ${MANIPULATOR_LIST}"Zpoly",\ )
    ADD_SUBDIRECTORY(Zpoly)
ENDIF(BUILD_PYRANHA_ZPOLY)

IF(BUILD_PYRANHA_QPOLY)
    SET(MODULE_LIST ${MODULE_LIST}"Qpoly",\ )
    SET(MANIPULATOR_LIST ${MANIPULATOR_LIST}"Qpoly",\ )
    ADD_SUBDIRECTORY(Qpoly)
ENDIF(BUILD_PYRANHA_QPOLY)

IF(BUILD_PYRANHA_QQPOLY)
    SET(MODULE_LIST ${MODULE_LIST}"Qqpoly",\ )
    SET(MANIPULATOR_LIST ${MANIPULATOR_LIST}"Qqpoly",\ )
    ADD_SUBDIRECTORY(Qqpoly)
ENDIF(BUILD_PYRANHA_QQPOLY)

IF(BUILD_PYRANHA_DPS)
    SET(MODULE_LIST ${MODULE_LIST}"Dps",\ )
    SET(MANIPULATOR_LIST ${MANIPULATOR_LIST}"Dps",\ )
    ADD_SUBDIRECTORY(Dps)
ENDIF(BUILD_PYRANHA_DPS)

IF(BUILD_PYRANHA_QPS)
    SET(MODULE_LIST ${MODULE_LIST}"Qps",\ )
    SET(MANIPULATOR_LIST ${MANIPULATOR_LIST}"Qps",\ )
    ADD_SUBDIRECTORY(Qps)
ENDIF(BUILD_PYRANHA_QPS)

IF(BUILD_PYRANHA_DQPS)
    SET(MODULE_LIST ${MODULE_LIST}"Dqps",\ )
    SET(MANIPULATOR_LIST ${MANIPULATOR_LIST}"Dqps",\ )
    ADD_SUBDIRECTORY(Dqps)
ENDIF(BUILD_PYRANHA_DQPS)

IF(BUILD_PYRANHA_QQPS)
    SET(MODULE_LIST ${MODULE_LIST}"Qqps",\ )
    SET(MANIPULATOR_LIST ${MANIPULATOR_LIST}"Qqps",\ )
    ADD_SUBDIRECTORY(Qqps)
ENDIF(BUILD_PYRANHA_QQPS)

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/template.global.__init__.py  ${CMAKE_CURRENT_SOURCE_DIR}/global.__init__.py)

INSTALL(FILES detail.py DESTINATION ${PYRANHA_INSTALL_PATH})
INSTALL(FILES global.__init__.py RENAME __init__.py DESTINATION ${PYRANHA_INSTALL_PATH})
