PROJECT(Piranha)

# Modified to only support MSVC >=2017 any longer and nothing else.
# The CMAKE is kept to keep some of the generations/installations easier
CMAKE_MINIMUM_REQUIRED(VERSION 2.6.0)

# This syntax was introduced by CMake 2.6 to handle backwards compatibility.
IF(COMMAND cmake_policy)
	CMAKE_POLICY(VERSION 2.6)
ENDIF(COMMAND cmake_policy)

ENABLE_LANGUAGE(CXX)
ENABLE_TESTING()


SET(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules)
INCLUDE(CMakeDependentOption)
INCLUDE(PiranhaCompilerLinkerSettings)
#INCLUDE(PiranhaPythonSetup)

INCLUDE(PiranhaVersionSetup)
# Place this after having determined the version number.
INCLUDE(PiranhaCPackSetup)


INCLUDE_DIRECTORIES("external_headers") # do we need them if we replace Boost with an external version??

# Main build options
OPTION(BUILD_ENABLE_ASSERTS "Enable asserts in code, for debugging." ON)

OPTION(BUILD_PYRANHA "Build Python bindings." OFF)
OPTION(BUILD_PIRANHA_STATIC "Build static Piranha library." OFF)
CMAKE_DEPENDENT_OPTION(BUILD_PYRANHA_DFS "Build double coefficient Fourier series manipulator in Pyranha." ON "BUILD_PYRANHA" OFF)
CMAKE_DEPENDENT_OPTION(BUILD_PYRANHA_DPOLY "Build double coefficient multivariate polynomial manipulator in Pyranha." ON "BUILD_PYRANHA" OFF)
CMAKE_DEPENDENT_OPTION(BUILD_PYRANHA_ZPOLY "Build arbitrary-size integer coefficient multivariate polynomial manipulator in Pyranha." ON "BUILD_PYRANHA" OFF)
CMAKE_DEPENDENT_OPTION(BUILD_PYRANHA_QPOLY "Build arbitrary-size rational coefficient multivariate polynomial manipulator in Pyranha." ON "BUILD_PYRANHA" OFF)
CMAKE_DEPENDENT_OPTION(BUILD_PYRANHA_QQPOLY "Build arbitrary-size rational coefficient multivariate polynomial manipulator with rational exponents in Pyranha." ON "BUILD_PYRANHA" OFF)
CMAKE_DEPENDENT_OPTION(BUILD_PYRANHA_DPS "Build double precision Poisson series manipulator in Pyranha." ON "BUILD_PYRANHA" OFF)
CMAKE_DEPENDENT_OPTION(BUILD_PYRANHA_QPS "Build rational coeffcient Poisson series manipulator in Pyranha." ON "BUILD_PYRANHA" OFF)
CMAKE_DEPENDENT_OPTION(BUILD_PYRANHA_DQPS "Build double coefficient - rational exponents Poisson series manipulator in Pyranha." ON "BUILD_PYRANHA" OFF)
CMAKE_DEPENDENT_OPTION(BUILD_PYRANHA_QQPS "Build rational coefficient - rational exponents Poisson series manipulator in Pyranha." ON "BUILD_PYRANHA" OFF)
CMAKE_DEPENDENT_OPTION(BUILD_PYRANHA_GUI "Build PyQt4 gui for pyranha." OFF "BUILD_PYRANHA" OFF) # TODO:?? upgrade to qt5??

#this also builds the static library. Should we have a specific option for the static library?
OPTION(BUILD_TESTS "Build performance and correctness tests." OFF)


# control build of load library (dll)        
IF(BUILD_PYRANHA)
        SET(BUILD_PIRANHA_DLL TRUE)
        INCLUDE(PiranhaPythonSetup)  # we need the python includes for piranha.dll
        PIRANHA_PYTHON_SETUP()
ENDIF()  

IF(BUILD_TESTS)
        SET(BUILD_PIRANHA_STATIC TRUE)
ENDIF()      
   
# Set default build type to "Release", change it in the GUI if you need to build with debug.
IF(NOT CMAKE_BUILD_TYPE)
	SET(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

SET(MANDATORY_LIBRARIES "")
SET(PIRANHA_DEFINITIONS "")

# Check for 64 bit support via pointer size.
IF(CMAKE_SIZEOF_VOID_P EQUAL 8)
        MESSAGE(STATUS "64 bit support detected")
        ADD_DEFINITIONS(-D_PIRANHA_64BIT)
ELSE()
        MESSAGE(FATAL_ERROR "No 64 bit support detected")
ENDIF()

#Cache size.
SET(CACHE_SIZE 512 CACHE STRING "Cache size in kilobytes.")
SET(PIRANHA_DEFINITIONS ${PIRANHA_DEFINITIONS} -D_PIRANHA_CACHE_SIZE=${CACHE_SIZE})
MARK_AS_ADVANCED(CACHE_SIZE)

MESSAGE(STATUS "Source directory: ${CMAKE_SOURCE_DIR}")

# Platform switches.
SET(PIRANHA_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})

IF(MSVC) 
        IF(MSVC_VERSION LESS 1910)
             MESSAGE(STATUS "This version ${MSVC_VERSION} of MSVC is not tested. It may work or not")
        ENDIF()
        SET(PIRANHA_DEFINITIONS ${PIRANHA_DEFINITIONS} -D__GMP_LIBGMP_DLL)                # i.e. it has to be the dll of MPIR
        
        SET(LIB_INSTALL_PATH .)
        SET(PIRANHA_INSTALL_PREFIX "${PIRANHA_INSTALL_PREFIX} ${PIRANHA_VERSION}")
ELSE(MSVC)
        MESSAGE(FATAL_ERROR "Only Microsoft Visual Studio >= 15 (2017) is supported")
ENDIF(MSVC)

# Find MPIR, the Windows incarnation of GMP/MPFR. The C and CXX are combined into the same dll and use the same header
FIND_PACKAGE(MPIR MODULE REQUIRED) # we have findMPIR in our local cmake script folder

# Note: the order seems important here, at least in MinGW.
SET(MANDATORY_LIBRARIES ${MANDATORY_LIBRARIES} ${MPIR_LIBRARIES})
INCLUDE_DIRECTORIES(${MPIR_INCLUDE_DIR})

# Common stuff
SET(THEORIES_INSTALL_PATH "theories")
SET(EXAMPLES_INSTALL_PATH "examples")

IF(BUILD_ENABLE_ASSERTS)
        SET(PIRANHA_DEFINITIONS ${PIRANHA_DEFINITIONS} -D_PIRANHA_ENABLE_ASSERTS)
ENDIF()

ADD_DEFINITIONS(${PIRANHA_DEFINITIONS})

# Add library directory
ADD_SUBDIRECTORY(src)

IF(BUILD_TESTS)
        ADD_SUBDIRECTORY(tests)
ENDIF(BUILD_TESTS)

# should be set form outside
SET(MY_PIRANHA_INSTALL "D:\\temp for tests\\Piranha")
# Installation directory
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
        SET(CMAKE_INSTALL_PREFIX ${MY_PIRANHA_INSTALL} CACHE PATH "..." FORCE)
endif()

# Pyranha: piranha's python bindings
IF(BUILD_PYRANHA)
            #$$	IF(WIN32)
            #$$		SET(WIN_PKG_ROOT ${CMAKE_SOURCE_DIR}/windows_package)
            #$$		SET(WIN_PKG_FILES
            #$$			${WIN_PKG_ROOT}/Console.exe
            #$$			${WIN_PKG_ROOT}/ConsoleHook.dll
            #$$			${WIN_PKG_ROOT}/FreeImage.dll
            #$$			${WIN_PKG_ROOT}/FreeImagePlus.dll
            #$$			${WIN_PKG_ROOT}/console.chm
            #$$			${WIN_PKG_ROOT}/console.xml
            #$$			${WIN_PKG_ROOT}/msvcp71.dll
            #$$			${WIN_PKG_ROOT}/msvcr71.dll
            #$$#			${CMAKE_SOURCE_DIR}/extras/ipy_profile_pyranha.py
            #$$#			${CMAKE_SOURCE_DIR}/extras/matplotlibrc
            #$$		)
            #$$		IF(MINGW)
            #$$			SET(WIN_PKG_FILES ${WIN_PKG_FILES} ${WIN_PKG_ROOT}/libgmp-3.dll ${WIN_PKG_ROOT}/libgmpxx-4.dll)
            #$$		ELSE(MINGW)
            #$$			SET(WIN_PKG_FILES ${WIN_PKG_FILES} ${WIN_PKG_ROOT}/mpir.dll)
            #$$		ENDIF(MINGW)
            #$$		IF(NOT MINGW)
            #$$			SET(WIN_MSVC_DLL_PATH ${CMAKE_SOURCE_DIR}/windows_package/Microsoft.VC90.CRT)
            #$$			SET(WIN_MSVC_DLL_FILES
            #$$				${WIN_MSVC_DLL_PATH}/Microsoft.VC90.CRT.manifest
            #$$				${WIN_MSVC_DLL_PATH}/msvcm90.dll
            #$$				${WIN_MSVC_DLL_PATH}/msvcp90.dll
            #$$				${WIN_MSVC_DLL_PATH}/msvcr90.dll
            #$$			)
            #$$			INSTALL(FILES ${WIN_MSVC_DLL_FILES} DESTINATION Microsoft.VC90.CRT)
            #$$		ENDIF(NOT MINGW)
            #$$		INSTALL(FILES ${WIN_PKG_FILES} DESTINATION .)
            #$$		# Now let's setup npp.
            #$$		INSTALL(FILES
            #$$			${WIN_PKG_ROOT}/npp/config.model.xml
            #$$			${WIN_PKG_ROOT}/npp/notepad++.exe
            #$$			${WIN_PKG_ROOT}/npp/doLocalConf.xml
            #$$			${WIN_PKG_ROOT}/npp/shortcuts.xml
            #$$			${WIN_PKG_ROOT}/npp/stylers.model.xml
            #$$			${WIN_PKG_ROOT}/npp/langs.model.xml
            #$$			${WIN_PKG_ROOT}/npp/readme.txt
            #$$			${WIN_PKG_ROOT}/npp/change.log
            #$$			${WIN_PKG_ROOT}/npp/license.txt
            #$$			${WIN_PKG_ROOT}/npp/SciLexer.dll
            #$$			${WIN_PKG_ROOT}/npp/contextMenu.xml
            #$$			${WIN_PKG_ROOT}/npp/LINEDRAW.TTF
            #$$			DESTINATION npp
            #$$		)
            #$$		INSTALL(FILES
            #$$			${WIN_PKG_ROOT}/npp/plugins/FTP_synchronize.dll
            #$$			${WIN_PKG_ROOT}/npp/plugins/NppAutoIndent.dll
            #$$			${WIN_PKG_ROOT}/npp/plugins/NppDocShare.dll
            #$$			${WIN_PKG_ROOT}/npp/plugins/NppExport.dll
            #$$			${WIN_PKG_ROOT}/npp/plugins/docMonitor.dll
            #$$			${WIN_PKG_ROOT}/npp/plugins/mimeTools.unicode.dll
            #$$			DESTINATION npp/plugins
            #$$		)
            #$$		INSTALL(FILES
            #$$			${WIN_PKG_ROOT}/npp/plugins/doc/FTP_synchonize.ReadMe.txt
            #$$			DESTINATION npp/plugins/doc
            #$$		)
            #$$		INSTALL(FILES
            #$$			${WIN_PKG_ROOT}/npp/plugins/Config/NppAutoIndent.ini
            #$$			DESTINATION npp/plugins/Config
            #$$		)
            #$$		INSTALL(FILES
            #$$			${WIN_PKG_ROOT}/npp/plugins/APIs/actionscript.xml
            #$$			${WIN_PKG_ROOT}/npp/plugins/APIs/cpp.xml
            #$$			${WIN_PKG_ROOT}/npp/plugins/APIs/css.xml
            #$$			${WIN_PKG_ROOT}/npp/plugins/APIs/html.xml
            #$$			${WIN_PKG_ROOT}/npp/plugins/APIs/javascript.xml
            #$$			${WIN_PKG_ROOT}/npp/plugins/APIs/nsis.xml
            #$$			${WIN_PKG_ROOT}/npp/plugins/APIs/perl.xml
            #$$			${WIN_PKG_ROOT}/npp/plugins/APIs/python.xml
            #$$			${WIN_PKG_ROOT}/npp/plugins/APIs/sql.xml
            #$$			${WIN_PKG_ROOT}/npp/plugins/APIs/vb.xml
            #$$			${WIN_PKG_ROOT}/npp/plugins/APIs/xml.xml
            #$$			${WIN_PKG_ROOT}/npp/plugins/APIs/c.xml
            #$$			${WIN_PKG_ROOT}/npp/plugins/APIs/cs.xml
            #$$			${WIN_PKG_ROOT}/npp/plugins/APIs/fortran.xml
            #$$			${WIN_PKG_ROOT}/npp/plugins/APIs/java.xml
            #$$			${WIN_PKG_ROOT}/npp/plugins/APIs/lisp.xml
            #$$			${WIN_PKG_ROOT}/npp/plugins/APIs/pascal.xml
            #$$			${WIN_PKG_ROOT}/npp/plugins/APIs/php.xml
            #$$			${WIN_PKG_ROOT}/npp/plugins/APIs/rc.xml
            #$$			${WIN_PKG_ROOT}/npp/plugins/APIs/tex.xml
            #$$			${WIN_PKG_ROOT}/npp/plugins/APIs/vhdl.xml
            #$$			DESTINATION npp/plugins/APIs
            #$$		)
            #$$		
            #$$		INSTALL(FILES 
            #$$			${CMAKE_SOURCE_DIR}/extras/ipy_profile_pyranha.py
            #$$#			${CMAKE_SOURCE_DIR}/extras/matplotlibrc
            #$$			DESTINATION .ipython/profile_pyranha/startup
            #$$		)
            #$$		
            #$$		INSTALL(FILES 
            #$$			${CMAKE_SOURCE_DIR}/extras/matplotlibrc
            #$$			DESTINATION .
            #$$		)
            #$$
            #$$	
            #$$	ENDIF(WIN32)
            
# Mandatory linking.
        SET(MANDATORY_LIBRARIES ${MANDATORY_LIBRARIES} piranha)
        ADD_SUBDIRECTORY(pyranha)
ENDIF(BUILD_PYRANHA)

            #$$# Install theories of motion and examples.
            #$$SET(THEORIES_FILES
            #$$	theories_of_motion/elp1.dps
            #$$	theories_of_motion/elp2.dfs
            #$$	theories_of_motion/elp3.dfs
            #$$)
            #$$
            #$$SET(EXAMPLES_FILES
            #$$	pyranha/examples/alt_besselJ.py
            #$$)
            #$$
            #$$INSTALL(FILES ${THEORIES_FILES} DESTINATION "${THEORIES_INSTALL_PATH}")
            #$$INSTALL(FILES ${EXAMPLES_FILES} DESTINATION "${EXAMPLES_INSTALL_PATH}")
            #$$IF(WIN32)
            #$$	CONFIGURE_FILE("${CMAKE_SOURCE_DIR}/COPYING" "${CMAKE_SOURCE_DIR}/windows_package/license.txt")
            #$$	CONFIGURE_FILE("${CMAKE_SOURCE_DIR}/Changelog" "${CMAKE_SOURCE_DIR}/windows_package/changelog.txt")
            #$$	INSTALL(FILES "${CMAKE_SOURCE_DIR}/windows_package/license.txt" "${CMAKE_SOURCE_DIR}/windows_package/changelog.txt" DESTINATION .)
            #$$ENDIF(WIN32)
            #$$
            #$$INCLUDE(CPack)
            #$$
            #$$