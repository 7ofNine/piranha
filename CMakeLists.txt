PROJECT(piranha)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6.0)

# This syntax was introduced by CMake 2.6 to handle backwards compatibility.
IF(COMMAND cmake_policy)
	CMAKE_POLICY(VERSION 2.6)
ENDIF(COMMAND cmake_policy)

SET(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules)
INCLUDE(CheckIncludeFiles)
INCLUDE(CMakeDependentOption)
INCLUDE(FindGMP)
INCLUDE(FindPythonLibs)
INCLUDE(PiranhaCPackSetup)
INCLUDE(PiranhaCompilerLinkerSettings)
INCLUDE(PiranhaPythonSetup)
INCLUDE(PiranhaVersionSetup)
ENABLE_LANGUAGE(CXX)
INCLUDE_DIRECTORIES("external_headers")

# Main build options
OPTION(BUILD_MAIN "Build 'main.cpp'." OFF)
OPTION(ENABLE_ASSERTS "Enable asserts in code, for debugging." ON)
OPTION(BUILD_PYRANHA "Build Python bindings." OFF)
CMAKE_DEPENDENT_OPTION(BUILD_PYRANHA_DFS "Build double coefficient Fourier series manipualtor in Pyranha." ON "BUILD_PYRANHA" OFF)
CMAKE_DEPENDENT_OPTION(BUILD_PYRANHA_DPOLY "Build double coefficient multivariate polynomial manipualtor in Pyranha." ON "BUILD_PYRANHA" OFF)
CMAKE_DEPENDENT_OPTION(BUILD_PYRANHA_ZPOLY "Build arbitrary-size integer coefficient multivariate polynomial manipualtor in Pyranha." ON "BUILD_PYRANHA" OFF)
CMAKE_DEPENDENT_OPTION(BUILD_PYRANHA_QPOLY "Build arbitrary-size rational coefficient multivariate polynomial manipualtor in Pyranha." ON "BUILD_PYRANHA" OFF)
CMAKE_DEPENDENT_OPTION(BUILD_PYRANHA_DPS "Build double precision Poisson series manipualtor in Pyranha." ON "BUILD_PYRANHA" OFF)
CMAKE_DEPENDENT_OPTION(BUILD_PYRANHA_QPS "Build rational coeffcient Poisson series manipualtor in Pyranha." ON "BUILD_PYRANHA" OFF)
CMAKE_DEPENDENT_OPTION(BUILD_PYRANHA_GUI "Build PyQt4 gui for pyranha." OFF "BUILD_PYRANHA" OFF)
OPTION(BUILD_TBB_MULTITHREADING "Use Intel's Thread Building Blocks for multi-threaded parallelism." OFF)
OPTION(BUILD_TESTS "Build performance and correctness tests." OFF)
OPTION(ENABLE_SSE2 "Enable usage of SSE2 instructions (requires a SSE2-enabled CPU) - EXPERIMENTAL." OFF)

SET(MANDATORY_LIBRARIES "")
SET(PIRANHA_DEFINITIONS "")

IF(BUILD_TBB_MULTITHREADING)
	IF(WIN32)
		# This is irrelevant at the moment, since TBB does not compile yet in MinGW.
		SET(PIRANHA_DEFINITIONS ${PIRANHA_DEFINITIONS} -DUSE_WINTHREAD)
	ELSE(WIN32)
		SET(PIRANHA_DEFINITIONS ${PIRANHA_DEFINITIONS} -DUSE_PTHREAD)
		SET(MANDATORY_LIBRARIES ${MANDATORY_LIBRARIES} pthread dl)
	ENDIF(WIN32)
	SET(PIRANHA_DEFINITIONS ${PIRANHA_DEFINITIONS} -D_PIRANHA_TBB)
ENDIF(BUILD_TBB_MULTITHREADING)

PIRANHA_COMPILER_LINKER_SETTINGS()
IF(BUILD_PYRANHA)
	PIRANHA_PYTHON_SETUP()
ENDIF(BUILD_PYRANHA)

# Platform switches
IF(UNIX)
	SET(LIB_INSTALL_PATH "lib")
ELSE(UNIX)
	SET(PIRANHA_DEFINITIONS ${PIRANHA_DEFINITIONS} -D_PIRANHA_WIN32 -D__GMP_LIBGMP_DLL)
	SET(LIB_INSTALL_PATH .)
ENDIF(UNIX)

# Find GMP.
FIND_PACKAGE(GMP REQUIRED)
# Note: the order seems important here, at least in MinGW.
SET(MANDATORY_LIBRARIES ${MANDATORY_LIBRARIES} ${GMPXX_LIBRARIES} ${GMP_LIBRARIES})
INCLUDE_DIRECTORIES(${GMP_INCLUDE_DIR} ${GMPXX_INCLUDE_DIR})

# SSE2
IF(ENABLE_SSE2)
	SET(PIRANHA_DEFINITIONS ${PIRANHA_DEFINITIONS} -D_PIRANHA_SSE2)
ENDIF(ENABLE_SSE2)

# Common stuff
SET(CMAKE_BUILD_TYPE "Release")
# Find out motion theories files to be installed
IF(UNIX)
	SET(THEORIES_INSTALL_PATH "share/piranha/theories")
ELSE(UNIX)
	SET(THEORIES_INSTALL_PATH "theories")
ENDIF(UNIX)
IF(ENABLE_ASSERTS)
	SET(PIRANHA_DEFINITIONS ${PIRANHA_DEFINITIONS} -D_ENABLE_ASSERTS)
ENDIF(ENABLE_ASSERTS)
IF(BUILD_TBB_MULTITHREADING)
	SET(PIRANHA_DEFINITIONS ${PIRANHA_DEFINITIONS} -D_REENTRANT)
ENDIF(BUILD_TBB_MULTITHREADING)

# TBB is not supported ATM on MinGW.
IF(BUILD_TBB_MULTITHREADING AND MINGW)
	MESSAGE(FATAL_ERROR "Intel TBB cannot be compiled in MinGW at the present time.")
ENDIF(BUILD_TBB_MULTITHREADING AND MINGW)

ADD_DEFINITIONS(${PIRANHA_DEFINITIONS})

# Add library directory
ADD_SUBDIRECTORY(src)

# Build a main.cpp file included in the source dir
IF(BUILD_MAIN)
	ADD_EXECUTABLE(main main.cpp)
	TARGET_LINK_LIBRARIES(main ${MANDATORY_LIBRARIES} piranha_static)
ENDIF(BUILD_MAIN)

IF(BUILD_TESTS)
	ADD_SUBDIRECTORY(tests)
ENDIF(BUILD_TESTS)

# pyranha: piranha's python bindings
IF(BUILD_PYRANHA)
	# Mandatory linking.
	SET(MANDATORY_LIBRARIES ${MANDATORY_LIBRARIES} piranha)
	ADD_SUBDIRECTORY(pyranha)
	IF(WIN32)
		SET(WIN_PYRANHA_PKG_FILES
			package/Console.exe
			package/ConsoleHook.dll
			package/FreeImage.dll
			package/FreeImagePlus.dll
			package/console.chm
			package/console.xml
			package/gmp.dll
			package/msvcp71.dll
			package/msvcr71.dll
			extras/ipy_profile_pyranha.py
		)
		INSTALL(FILES ${WIN_PYRANHA_PKG_FILES} DESTINATION .)
	ENDIF(WIN32)
ENDIF(BUILD_PYRANHA)

# Install phase: find and install theories of motion
FILE(GLOB_RECURSE THEORIES theories_of_motion/*.csv theories_of_motion/*.phl theories_of_motion/*.dat)
INSTALL(FILES ${THEORIES} DESTINATION ${THEORIES_INSTALL_PATH})

INCLUDE(CPack)
