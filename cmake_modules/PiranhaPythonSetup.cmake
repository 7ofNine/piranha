MACRO(PIRANHA_PYTHON_SETUP)
  IF(BUILD_PYRANHA)
    # Find Python libraries
    INCLUDE(${CMAKE_ROOT}/Modules/FindPythonLibs.cmake)
    FIND_PACKAGE(PythonLibs REQUIRED)
    INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
    MESSAGE(STATUS "Python library; " "${PYTHON_LIBRARY}")
    STRING(REGEX MATCH libpython[0-9]\\.?[0-9] PYTHON_LIBRARY_VERSION ${PYTHON_LIBRARY})
    STRING(REGEX REPLACE \\. "" PYTHON_LIBRARY_VERSION ${PYTHON_LIBRARY_VERSION})
    STRING(REGEX REPLACE libpython "" PYTHON_LIBRARY_VERSION ${PYTHON_LIBRARY_VERSION})
    MESSAGE(STATUS "Python library version: " ${PYTHON_LIBRARY_VERSION})
    SET(PYDEXTENSION FALSE)
    IF(${PYTHON_LIBRARY_VERSION} GREATER 24 AND WIN32)
      MESSAGE(STATUS "Python >= 2.5 detected on Windows platform. Output extension for compiled modules will be '.pyd'.")
      SET(PYDEXTENSION TRUE)
    ENDIF(${PYTHON_LIBRARY_VERSION} GREATER 24 AND WIN32)
    # Trick to locate python's modules directory
    IF(WIN32)
      STRING(REGEX REPLACE "libs/libpython.*" "DLLs/" PYTHON_MODULES_PATH ${PYTHON_LIBRARIES})
    ELSE(WIN32)
      STRING(REGEX REPLACE "config/libpython.*" "site-packages/" PYTHON_MODULES_PATH ${PYTHON_LIBRARIES})
    ENDIF(WIN32)
    MESSAGE(STATUS "Python modules install path: " "${PYTHON_MODULES_PATH}")
  ENDIF(BUILD_PYRANHA)
ENDMACRO(PIRANHA_PYTHON_SETUP)