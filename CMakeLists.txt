#
# mainline/CMakeList
#
# Modified to only support MSVC >=2019 with C++20 any longer and nothing else.
# This is old style unstructured CMAKE!
#

CMAKE_MINIMUM_REQUIRED(VERSION 3.20.0 FATAL_ERROR)

PROJECT(Piranha CXX)

# we want to organize projects in folders
SET_PROPERTY(GLOBAL PROPERTY USE_FOLDERS ON)

set(CMAKE_CXX_STANDARD 20)     # we require C++20 even if we don't use many features, yet! But that may come one day

# we only support MSVC >=2019 and 64 Bit
IF(MSVC)
        IF(MSVC_VERSION LESS 1920)
             MESSAGE(STATUS "This version ${MSVC_VERSION} of MSVC is not tested. It may work or not")
        ENDIF()

        # Check for 64 bit support via pointer size.
        IF(CMAKE_SIZEOF_VOID_P EQUAL 8)
            MESSAGE(STATUS "64 bit support detected")
        ELSE()
            MESSAGE(FATAL_ERROR "No 64 bit support detected")
        ENDIF()

ELSE(MSVC)
        MESSAGE(FATAL_ERROR "Only Microsoft Visual Studio >= 16 (2019) is supported")
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

CMAKE_DEPENDENT_OPTION(BUILD_PYRANHA_DFS "Build double coefficient Fourier series manipulator in Pyranha." ON "BUILD_PYRANHA" OFF)
CMAKE_DEPENDENT_OPTION(BUILD_PYRANHA_DPOLY "Build double coefficient multivariate polynomial manipulator in Pyranha." ON "BUILD_PYRANHA" OFF)
CMAKE_DEPENDENT_OPTION(BUILD_PYRANHA_ZPOLY "Build arbitrary-size integer coefficient multivariate polynomial manipulator in Pyranha." ON "BUILD_PYRANHA" OFF)
CMAKE_DEPENDENT_OPTION(BUILD_PYRANHA_QPOLY "Build arbitrary-size rational coefficient multivariate polynomial manipulator in Pyranha." ON "BUILD_PYRANHA" OFF)
CMAKE_DEPENDENT_OPTION(BUILD_PYRANHA_QQPOLY "Build arbitrary-size rational coefficient multivariate polynomial manipulator with rational exponents in Pyranha." ON "BUILD_PYRANHA" OFF)
CMAKE_DEPENDENT_OPTION(BUILD_PYRANHA_DPS "Build double precision Poisson series manipulator in Pyranha." ON "BUILD_PYRANHA" OFF)
CMAKE_DEPENDENT_OPTION(BUILD_PYRANHA_QPS "Build rational coeffcient Poisson series manipulator in Pyranha." ON "BUILD_PYRANHA" OFF)
CMAKE_DEPENDENT_OPTION(BUILD_PYRANHA_DQPS "Build double coefficient - rational exponents Poisson series manipulator in Pyranha." ON "BUILD_PYRANHA" OFF)
CMAKE_DEPENDENT_OPTION(BUILD_PYRANHA_QQPS "Build rational coefficient - rational exponents Poisson series manipulator in Pyranha." ON "BUILD_PYRANHA" OFF)
CMAKE_DEPENDENT_OPTION(BUILD_PYRANHA_GUI "Build PySide2 gui for pyranha." OFF "BUILD_PYRANHA" OFF) # TODO:?? upgrade to qt5??

OPTION(BUILD_TESTS "Build performance and correctness tests." OFF)




# Python setup
IF(BUILD_PYRANHA) 
    INCLUDE(PiranhaPythonSetup) # we need the python includes for piranha.dll ?
    PIRANHA_PYTHON_SETUP()
ENDIF()

#Boost setup
INCLUDE(PiranhaFindBoost)
PIRANHA_BOOST_SETUP()

MESSAGE(STATUS "thread library: ${Boost_thread_LIBRARY_RELEASE}")

# Find MPIR, the Windows incarnation of GMP/MPFR. The C and CXX are combined into the same dll and use the same header
FIND_PACKAGE(MPIR MODULE REQUIRED) # we have findMPIR in our local cmake script folder




# Set default build type to "Release", change it in the VS GUI if you need to build with debug.
IF(NOT CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

SET(MANDATORY_LIBRARIES "")
SET(PIRANHA_DEFINITIONS "")



#Cache size.
SET(CACHE_SIZE 1024 CACHE STRING "Cache size in kilobytes.")
SET(PIRANHA_DEFINITIONS ${PIRANHA_DEFINITIONS} PIRANHA_CACHE_SIZE=${CACHE_SIZE})
MARK_AS_ADVANCED(CACHE_SIZE)

#Build with asserts
IF(BUILD_ENABLE_ASSERTS)
        SET(PIRANHA_DEFINITIONS ${PIRANHA_DEFINITIONS} PIRANHA_ENABLE_ASSERTS)
ENDIF()

ADD_COMPILE_DEFINITIONS(${PIRANHA_DEFINITIONS})


MESSAGE(STATUS "Source directory: ${CMAKE_SOURCE_DIR}")

# Platform switches.
SET(PIRANHA_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})



# Note: the order seems important here, at least in MinGW it was
SET(MANDATORY_LIBRARIES ${MANDATORY_LIBRARIES} ${MPIR_LIBRARIES})
INCLUDE_DIRECTORIES(${MPIR_INCLUDE_DIR}) # MPIR is probably needed evwerywhere



# Add libraries piranha.dll
ADD_SUBDIRECTORY(src)

# Add the tests
IF(BUILD_TESTS)
        ADD_SUBDIRECTORY(tests)
ENDIF(BUILD_TESTS)

SET(PIRANHA_INSTALL_PATH ${CMAKE_INSTALL_PREFIX} CACHE PATH "Path to the PYRANHA installation" FORCE)


# Pyranha: piranha's python bindings
IF(BUILD_PYRANHA)
    SET(MANDATORY_LIBRARIES ${MANDATORY_LIBRARIES} piranha) # the piranha.dll is needed
    ADD_SUBDIRECTORY(pyranha) #uses MANDATORY_LIBRARIES

# TEST
#TEST
SET(TEST_STRING $<TARGET_FILE:Boost::thread>)
MESSAGE(STATUS "Boost thread: ${TEST_STRING}")

#TEST
#TEST	
    
    
	# installation of PYRANHA
	SET(WIN_PKG_ROOT ${CMAKE_SOURCE_DIR}/WinPackage)
	MESSAGE (STATUS "winpackage: ${WIN_PKG_ROOT}") 
    
    SET(WIN_PKG_FILES
						${WIN_PKG_ROOT}/Console.exe
						${WIN_PKG_ROOT}/ConsoleHook.dll
						${WIN_PKG_ROOT}/FreeImagePlus.dll
            #$$         ${WIN_PKG_ROOT}/console.chm   # TODO: get it from somewhere or build it
                        ${WIN_PKG_ROOT}/license.txt
                        ${CMAKE_SOURCE_DIR}/extras/matplotlibrc
		)
	SET(WIN_PKG_FILES ${WIN_PKG_FILES} ${MPIR_INCLUDE_DIR}/mpir.dll) # header in library for mpir are in the same folder
    
    INSTALL(FILES ${WIN_PKG_FILES} DESTINATION .)  # DESTINATION is relative to CMAKE_INSTALL_PREFIX

    #files coming from a different source or have to be configured
    CONFIGURE_FILE(${WIN_PKG_ROOT}/template.console.xml ${WIN_PKG_ROOT}/generated.console.xml) #Contains the path to python scripts, i. e. has to be after Python setup 
	INSTALL(FILES ${WIN_PKG_ROOT}/generated.console.xml RENAME console.xml DESTINATION .) #copy configuration for console. 
    
    CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/extras/template.ipy_profile_pyranha.py ${CMAKE_CURRENT_SOURCE_DIR}/extras/ipy_profile_pyranha.py)
	INSTALL(FILES  ${CMAKE_SOURCE_DIR}/extras/ipy_profile_pyranha.py  DESTINATION .ipython/profile_pyranha/startup) # install pyranha startup for ipython

    CONFIGURE_FILE("${CMAKE_SOURCE_DIR}/COPYING" "${WIN_PKG_ROOT}/license.txt")   # in case we have to have something to make dependent on build
	INSTALL(FILES "${WIN_PKG_ROOT}/license.txt" DESTINATION .)
    
    INSTALL(FILES ${CMAKE_SOURCE_DIR}/extras/matplotlibrc DESTINATION .)  # why is it installed in two different locations???? something wrong with configuration
ENDIF(BUILD_PYRANHA)


# Common installation files
SET(THEORIES_INSTALL_PATH "theories")   # relative within destination see below
SET(EXAMPLES_INSTALL_PATH "examples")
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
            
#CONFIGURE_FILE("${CMAKE_SOURCE_DIR}/COPYING" "${WIN_PKG_ROOT}/license.txt")

#CONFIGURE_FILE("${CMAKE_SOURCE_DIR}/Changelog" "${CMAKE_SOURCE_DIR}/windows_package/changelog.txt")   #ond day we update this one
#INSTALL(FILES "${CMAKE_SOURCE_DIR}/windows_package/license.txt" "${CMAKE_SOURCE_DIR}/windows_package/changelog.txt" DESTINATION .) # change to winPackage

#$$INCLUDE(CPack)
