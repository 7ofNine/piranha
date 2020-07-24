# Modified to only support MSVC >=2017 any longer and nothing else.
# The CMAKE is kept to keep some of the generations/installations easier
#
# This is old style unstructured CMAKE!
#
CMAKE_MINIMUM_REQUIRED(VERSION 3.10.0 FATAL_ERROR)

PROJECT(Piranha CXX)

# This syntax was introduced by CMake 2.6 to handle backwards compatibility.
IF(COMMAND cmake_policy)
    CMAKE_POLICY(VERSION 2.6)
ENDIF(COMMAND cmake_policy)

# we want to organize projects in folders
SET_PROPERTY(GLOBAL PROPERTY USE_FOLDERS ON)

# we only support MSVC >=2017 and 64 Bit
IF(MSVC)
        set(CMAKE_CXX_STANDARD 17) # 17 still has warnigns and errors
        IF(MSVC_VERSION LESS 1910)
             MESSAGE(STATUS "This version ${MSVC_VERSION} of MSVC is not tested. It may work or not")
        ENDIF()

        # Check for 64 bit support via pointer size.
        IF(CMAKE_SIZEOF_VOID_P EQUAL 8)
            MESSAGE(STATUS "64 bit support detected")
        ELSE()
            MESSAGE(FATAL_ERROR "No 64 bit support detected")
        ENDIF()

#        SET(PIRANHA_DEFINITIONS ${PIRANHA_DEFINITIONS} -D__GMP_LIBGMP_DLL)                # i.e. it has to be the dll of MPIR, do we really need that?

#        SET(LIB_INSTALL_PATH .)
#        SET(PIRANHA_INSTALL_PREFIX "${PIRANHA_INSTALL_PREFIX} ${PIRANHA_VERSION}")
ELSE(MSVC)
        MESSAGE(FATAL_ERROR "Only Microsoft Visual Studio >= 15 (2017) is supported")
ENDIF(MSVC)


ENABLE_TESTING()

SET(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules) # Our own CMAKE scripts
INCLUDE(CMakeDependentOption)
INCLUDE(PiranhaCompilerLinkerSettings)
INCLUDE(PiranhaVersionSetup)

# Place this after having determined the version number.
INCLUDE(PiranhaCPackSetup)


# Main build options
OPTION(BUILD_ENABLE_ASSERTS "Enable asserts in code, for debugging." ON)

OPTION(BUILD_PYRANHA "Build Python bindings." OFF)
#PYTHON 2 no longer supported
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

OPTION(BUILD_TESTS "Build performance and correctness tests." OFF)




# control build of load library (dll)
IF(BUILD_PYRANHA) # ultimately that should always be 3.* . Temporary to be able to switch between old and new
    SET(PYRANHA_PYTHON_BUILD_VERSION 3)
    SET(PYRANHA_BOOST_PYTHON_TARGET "python38") # depdends on the precise version of python should be set up depdedent on python)
    INCLUDE(PiranhaPythonSetup) # we need the python includes for piranha.dll ?
    PIRANHA_PYTHON_SETUP(${PYRANHA_PYTHON_BUILD_VERSION}) 
ENDIF()

# Set default build type to "Release", change it in the GUI if you need to build with debug.
IF(NOT CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

SET(MANDATORY_LIBRARIES "")
SET(PIRANHA_DEFINITIONS "")


#Cache size.
SET(CACHE_SIZE 1024 CACHE STRING "Cache size in kilobytes.")
SET(PIRANHA_DEFINITIONS ${PIRANHA_DEFINITIONS} -D_PIRANHA_CACHE_SIZE=${CACHE_SIZE})
MARK_AS_ADVANCED(CACHE_SIZE)

IF(BUILD_ENABLE_ASSERTS)
        SET(PIRANHA_DEFINITIONS ${PIRANHA_DEFINITIONS} -D_PIRANHA_ENABLE_ASSERTS)
ENDIF()

ADD_DEFINITIONS(${PIRANHA_DEFINITIONS})





MESSAGE(STATUS "Source directory: ${CMAKE_SOURCE_DIR}")

# Platform switches.
SET(PIRANHA_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})


# Find MPIR, the Windows incarnation of GMP/MPFR. The C and CXX are combined into the same dll and use the same header
FIND_PACKAGE(MPIR MODULE REQUIRED) # we have findMPIR in our local cmake script folder


# Note: the order seems important here, at least in MinGW.
SET(MANDATORY_LIBRARIES ${MANDATORY_LIBRARIES} ${MPIR_LIBRARIES})
INCLUDE_DIRECTORIES(${MPIR_INCLUDE_DIR}) # MPIR is probably needed evwerywhere

# Common stuff
SET(THEORIES_INSTALL_PATH "theories")   # relative within destination see below
SET(EXAMPLES_INSTALL_PATH "examples")


# Add libraries DLL
ADD_SUBDIRECTORY(src)

# Add the tests
IF(BUILD_TESTS)
        ADD_SUBDIRECTORY(tests)
ENDIF(BUILD_TESTS)

# should be set from outside
SET(MY_PIRANHA_INSTALL "D:\\temp for tests\\${CMAKE_PROJECT_NAME}test")  # temporarily hardcoded for testing purposes
# Installation directory
#if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
        SET(CMAKE_INSTALL_PREFIX ${MY_PIRANHA_INSTALL} CACHE PATH "Path to the PYRANHA installation" FORCE)
#endif()

# Pyranha: piranha's python bindings
IF(BUILD_PYRANHA)
    SET(MANDATORY_LIBRARIES ${MANDATORY_LIBRARIES} piranha) # the piranha.dll is needed
    ADD_SUBDIRECTORY(pyranha)
	
	# installation of PYRANHA
	SET(WIN_PKG_ROOT ${CMAKE_SOURCE_DIR}/WinPackage)
#	CONFIGURE_FILE(${WIN_PKG_ROOT}/template.console.xml ${WIN_PKG_ROOT}/generated.console.xml)  #TODO: one day we have to change to consoleZ
	MESSAGE (STATUS "winpackage: ${WIN_PKG_ROOT}") 
    SET(WIN_PKG_FILES
						${WIN_PKG_ROOT}/Console.exe
						${WIN_PKG_ROOT}/ConsoleHook.dll
						${WIN_PKG_ROOT}/FreeImagePlus.dll
            #$$         ${WIN_PKG_ROOT}/console.chm   # TODO: get it from somewhere or build it
						${WIN_PKG_ROOT}/console.xml
            #            ${CMAKE_SOURCE_DIR}/extras/ipy_profile_pyranha.py
                        ${CMAKE_SOURCE_DIR}/extras/matplotlibrc
		)
	MESSAGE (STATUS "${WIN_PKG_ROOT}/console.xml")	
	SET(WIN_PKG_FILES ${WIN_PKG_FILES} ${MPIR_INCLUDE_DIR}/mpir.dll) # header in library for mpir are in the same folder
		
    INSTALL(FILES ${WIN_PKG_FILES} DESTINATION .)  # DESTINATION is relative to CMAKE_INSTALL_PREFIX
	
	INSTALL(FILES  ${CMAKE_SOURCE_DIR}/extras/ipy_profile_pyranha.py  DESTINATION .ipython/profile_pyranha/startup) # install pyranha startup for ipython

	
	## TODO: npp. we have it installed. should find it and not reinstall whatever outdated version of it
            #$$     # Now let's setup npp. TODO: do we need that?
            #$$     INSTALL(FILES
            #$$         ${WIN_PKG_ROOT}/npp/config.model.xml
            #$$         ${WIN_PKG_ROOT}/npp/notepad++.exe
            #$$         ${WIN_PKG_ROOT}/npp/doLocalConf.xml
            #$$         ${WIN_PKG_ROOT}/npp/shortcuts.xml
            #$$         ${WIN_PKG_ROOT}/npp/stylers.model.xml
            #$$         ${WIN_PKG_ROOT}/npp/langs.model.xml
            #$$         ${WIN_PKG_ROOT}/npp/readme.txt
            #$$         ${WIN_PKG_ROOT}/npp/change.log
            #$$         ${WIN_PKG_ROOT}/npp/license.txt
            #$$         ${WIN_PKG_ROOT}/npp/SciLexer.dll
            #$$         ${WIN_PKG_ROOT}/npp/contextMenu.xml
            #$$         ${WIN_PKG_ROOT}/npp/LINEDRAW.TTF
            #$$         DESTINATION npp
            #$$     )
            #$$     INSTALL(FILES
            #$$         ${WIN_PKG_ROOT}/npp/plugins/FTP_synchronize.dll
            #$$         ${WIN_PKG_ROOT}/npp/plugins/NppAutoIndent.dll
            #$$         ${WIN_PKG_ROOT}/npp/plugins/NppDocShare.dll
            #$$         ${WIN_PKG_ROOT}/npp/plugins/NppExport.dll
            #$$         ${WIN_PKG_ROOT}/npp/plugins/docMonitor.dll
            #$$         ${WIN_PKG_ROOT}/npp/plugins/mimeTools.unicode.dll
            #$$         DESTINATION npp/plugins
            #$$     )
            #$$     INSTALL(FILES
            #$$         ${WIN_PKG_ROOT}/npp/plugins/doc/FTP_synchonize.ReadMe.txt
            #$$         DESTINATION npp/plugins/doc
            #$$     )
            #$$     INSTALL(FILES
            #$$         ${WIN_PKG_ROOT}/npp/plugins/Config/NppAutoIndent.ini
            #$$         DESTINATION npp/plugins/Config
            #$$     )
            #$$     INSTALL(FILES
            #$$         ${WIN_PKG_ROOT}/npp/plugins/APIs/actionscript.xml
            #$$         ${WIN_PKG_ROOT}/npp/plugins/APIs/cpp.xml
            #$$         ${WIN_PKG_ROOT}/npp/plugins/APIs/css.xml
            #$$         ${WIN_PKG_ROOT}/npp/plugins/APIs/html.xml
            #$$         ${WIN_PKG_ROOT}/npp/plugins/APIs/javascript.xml
            #$$         ${WIN_PKG_ROOT}/npp/plugins/APIs/nsis.xml
            #$$         ${WIN_PKG_ROOT}/npp/plugins/APIs/perl.xml
            #$$         ${WIN_PKG_ROOT}/npp/plugins/APIs/python.xml
            #$$         ${WIN_PKG_ROOT}/npp/plugins/APIs/sql.xml
            #$$         ${WIN_PKG_ROOT}/npp/plugins/APIs/vb.xml
            #$$         ${WIN_PKG_ROOT}/npp/plugins/APIs/xml.xml
            #$$         ${WIN_PKG_ROOT}/npp/plugins/APIs/c.xml
            #$$         ${WIN_PKG_ROOT}/npp/plugins/APIs/cs.xml
            #$$         ${WIN_PKG_ROOT}/npp/plugins/APIs/fortran.xml
            #$$         ${WIN_PKG_ROOT}/npp/plugins/APIs/java.xml
            #$$         ${WIN_PKG_ROOT}/npp/plugins/APIs/lisp.xml
            #$$         ${WIN_PKG_ROOT}/npp/plugins/APIs/pascal.xml
            #$$         ${WIN_PKG_ROOT}/npp/plugins/APIs/php.xml
            #$$         ${WIN_PKG_ROOT}/npp/plugins/APIs/rc.xml
            #$$         ${WIN_PKG_ROOT}/npp/plugins/APIs/tex.xml
            #$$         ${WIN_PKG_ROOT}/npp/plugins/APIs/vhdl.xml
            #$$         DESTINATION npp/plugins/APIs
            #$$     )
            #$$
                 INSTALL(FILES
            #$$         ${CMAKE_SOURCE_DIR}/extras/ipy_profile_pyranha.py  # we already did that. Whyis it here again
                        ${CMAKE_SOURCE_DIR}/extras/matplotlibrc
                        DESTINATION .ipython/profile_pyranha/startup
                 )
            
                 INSTALL(FILES
                     ${CMAKE_SOURCE_DIR}/extras/matplotlibrc   # why is it installed in two different locations???? something wrong with configuration
                     DESTINATION .
                 )
            
            
            ## ENDIF(WIN32) #delete?
ENDIF(BUILD_PYRANHA)

            # Install theories of motion and examples.
            SET(THEORIES_FILES
                theories_of_motion/elp1.dps   # should be relative to source directory
                theories_of_motion/elp2.dfs
                theories_of_motion/elp3.dfs
            )
            
            SET(EXAMPLES_FILES
                pyranha/examples/alt_besselJ.py
            )
            
            INSTALL(FILES ${THEORIES_FILES} DESTINATION "${THEORIES_INSTALL_PATH}")
            INSTALL(FILES ${EXAMPLES_FILES} DESTINATION "${EXAMPLES_INSTALL_PATH}")
            
            
            CONFIGURE_FILE("${CMAKE_SOURCE_DIR}/COPYING" "${CMAKE_SOURCE_DIR}/windows_package/license.txt")
            CONFIGURE_FILE("${CMAKE_SOURCE_DIR}/Changelog" "${CMAKE_SOURCE_DIR}/windows_package/changelog.txt")
 #           INSTALL(FILES "${CMAKE_SOURCE_DIR}/windows_package/license.txt" "${CMAKE_SOURCE_DIR}/windows_package/changelog.txt" DESTINATION .) # change to winPackae
            
            #$$
            #$$INCLUDE(CPack)
            #$$
            #$$
