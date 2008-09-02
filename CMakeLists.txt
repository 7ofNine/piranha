PROJECT(Piranha)

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
INCLUDE(PiranhaCompilerLinkerSettings)
INCLUDE(PiranhaPythonSetup)
INCLUDE(PiranhaVersionSetup)
# Place this after having determined the version number.
INCLUDE(PiranhaCPackSetup)
ENABLE_LANGUAGE(CXX)
ENABLE_TESTING()
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
OPTION(BUILD_MULTITHREADING "Enable multi-threaded parallelism." OFF)
OPTION(BUILD_TESTS "Build performance and correctness tests." OFF)
OPTION(ENABLE_SSE2 "Enable usage of SSE2 instructions (requires a SSE2-enabled CPU) - EXPERIMENTAL." OFF)

SET(MANDATORY_LIBRARIES "")
SET(PIRANHA_DEFINITIONS "")

IF(BUILD_MULTITHREADING)
	IF(WIN32)
		# This is irrelevant at the moment, since TBB does not compile yet in MinGW.
		SET(PIRANHA_DEFINITIONS ${PIRANHA_DEFINITIONS} -DUSE_WINTHREAD)
	ELSE(WIN32)
		SET(PIRANHA_DEFINITIONS ${PIRANHA_DEFINITIONS} -DUSE_PTHREAD)
		SET(MANDATORY_LIBRARIES ${MANDATORY_LIBRARIES} pthread dl)
	ENDIF(WIN32)
	SET(PIRANHA_DEFINITIONS ${PIRANHA_DEFINITIONS} -D_PIRANHA_MT)
ENDIF(BUILD_MULTITHREADING)

PIRANHA_COMPILER_LINKER_SETTINGS()
IF(BUILD_PYRANHA)
	PIRANHA_PYTHON_SETUP()
ENDIF(BUILD_PYRANHA)

# Platform switches
SET(PIRANHA_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
IF(UNIX)
	SET(LIB_INSTALL_PATH "lib")
ELSE(UNIX)
	SET(PIRANHA_DEFINITIONS ${PIRANHA_DEFINITIONS} -D_PIRANHA_WIN32)
	IF(MINGW)
		SET(PIRANHA_DEFINITIONS ${PIRANHA_DEFINITIONS} -D_PIRANHA_MINGW)
	ELSE(MINGW)
		SET(PIRANHA_DEFINITIONS ${PIRANHA_DEFINITIONS} -D__GMP_LIBGMP_DLL)
	ENDIF(MINGW)
	SET(LIB_INSTALL_PATH .)
	SET(PIRANHA_INSTALL_PREFIX "${PIRANHA_INSTALL_PREFIX} ${PIRANHA_VERSION}")
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
	SET(PIRANHA_DEFINITIONS ${PIRANHA_DEFINITIONS} -D_PIRANHA_ENABLE_ASSERTS)
ENDIF(ENABLE_ASSERTS)
IF(BUILD_MULTITHREADING)
	SET(PIRANHA_DEFINITIONS ${PIRANHA_DEFINITIONS} -D_REENTRANT)
ENDIF(BUILD_MULTITHREADING)

# TBB is not supported ATM on MinGW.
IF(BUILD_MULTITHREADING AND MINGW)
	MESSAGE(FATAL_ERROR "Intel's Thread Building Blocks cannot be compiled in MinGW at the present time.")
ENDIF(BUILD_MULTITHREADING AND MINGW)

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

# Pyranha: piranha's python bindings
IF(BUILD_PYRANHA)
	IF(WIN32)
		SET(WIN_PKG_ROOT ${CMAKE_SOURCE_DIR}/windows_package)
		SET(WIN_PKG_FILES
			${WIN_PKG_ROOT}/Console.exe
			${WIN_PKG_ROOT}/ConsoleHook.dll
			${WIN_PKG_ROOT}/FreeImage.dll
			${WIN_PKG_ROOT}/FreeImagePlus.dll
			${WIN_PKG_ROOT}/console.chm
			${WIN_PKG_ROOT}/console.xml
			${WIN_PKG_ROOT}/msvcp71.dll
			${WIN_PKG_ROOT}/msvcr71.dll
			${CMAKE_SOURCE_DIR}/extras/ipy_profile_pyranha.py
			${CMAKE_SOURCE_DIR}/extras/matplotlibrc
			${WIN_PKG_ROOT}/npp
			${WIN_PKG_ROOT}/npp/config.model.xml
			${WIN_PKG_ROOT}/npp/notepad++.exe
			${WIN_PKG_ROOT}/npp/doLocalConf.xml
			${WIN_PKG_ROOT}/npp/shortcuts.xml
			${WIN_PKG_ROOT}/npp/stylers.model.xml
			${WIN_PKG_ROOT}/npp/langs.model.xml
			${WIN_PKG_ROOT}/npp/readme.txt
			${WIN_PKG_ROOT}/npp/change.log
			${WIN_PKG_ROOT}/npp/license.txt
			${WIN_PKG_ROOT}/npp/SciLexer.dll
			${WIN_PKG_ROOT}/npp/contextMenu.xml
			${WIN_PKG_ROOT}/npp/plugins
			${WIN_PKG_ROOT}/npp/plugins/doc
			${WIN_PKG_ROOT}/npp/plugins/doc/NppExec.txt
			${WIN_PKG_ROOT}/npp/plugins/doc/NPPTextFXdemo.TXT
			${WIN_PKG_ROOT}/npp/plugins/doc/NppExec_TechInfo.txt
			${WIN_PKG_ROOT}/npp/plugins/doc/FTP_synchonize.ReadMe.txt
			${WIN_PKG_ROOT}/npp/plugins/APIs
			${WIN_PKG_ROOT}/npp/plugins/APIs/lisp.xml
			${WIN_PKG_ROOT}/npp/plugins/APIs/cpp.xml
			${WIN_PKG_ROOT}/npp/plugins/APIs/xml.xml
			${WIN_PKG_ROOT}/npp/plugins/APIs/css.xml
			${WIN_PKG_ROOT}/npp/plugins/APIs/c.xml
			${WIN_PKG_ROOT}/npp/plugins/APIs/actionscript.xml
			${WIN_PKG_ROOT}/npp/plugins/APIs/html.xml
			${WIN_PKG_ROOT}/npp/plugins/APIs/pascal.xml
			${WIN_PKG_ROOT}/npp/plugins/APIs/sql.xml
			${WIN_PKG_ROOT}/npp/plugins/APIs/tex.xml
			${WIN_PKG_ROOT}/npp/plugins/APIs/cs.xml
			${WIN_PKG_ROOT}/npp/plugins/APIs/python.xml
			${WIN_PKG_ROOT}/npp/plugins/APIs/nsis.xml
			${WIN_PKG_ROOT}/npp/plugins/APIs/rc.xml
			${WIN_PKG_ROOT}/npp/plugins/APIs/vb.xml
			${WIN_PKG_ROOT}/npp/plugins/APIs/java.xml
			${WIN_PKG_ROOT}/npp/plugins/APIs/vhdl.xml
			${WIN_PKG_ROOT}/npp/plugins/APIs/javascript.xml
			${WIN_PKG_ROOT}/npp/plugins/APIs/php.xml
			${WIN_PKG_ROOT}/npp/plugins/APIs/fortran.xml
			${WIN_PKG_ROOT}/npp/plugins/APIs/perl.xml
			${WIN_PKG_ROOT}/npp/plugins/Config
			${WIN_PKG_ROOT}/npp/plugins/mimeTools.dll
			${WIN_PKG_ROOT}/npp/plugins/docMonitor.dll
			${WIN_PKG_ROOT}/npp/plugins/NPPTextFX.dll
			${WIN_PKG_ROOT}/npp/plugins/NppExport.dll
			${WIN_PKG_ROOT}/npp/plugins/NppExec.dll
			${WIN_PKG_ROOT}/npp/plugins/LightExplorer.dll
			${WIN_PKG_ROOT}/npp/plugins/SpellChecker.dll
			${WIN_PKG_ROOT}/npp/plugins/ComparePlugin.dll
			${WIN_PKG_ROOT}/npp/plugins/FTP_synchronizeA.dll
			${WIN_PKG_ROOT}/npp/plugins/NPPTextFX
			${WIN_PKG_ROOT}/npp/plugins/NPPTextFX/libTidy.dll
			${WIN_PKG_ROOT}/npp/plugins/NPPTextFX/TIDYCFG.INI
			${WIN_PKG_ROOT}/npp/plugins/NPPTextFX/AsciiToEBCDIC.bin
			${WIN_PKG_ROOT}/npp/plugins/NPPTextFX/W3C-HTMLValidator.htm
			${WIN_PKG_ROOT}/npp/plugins/NPPTextFX/W3C-CSSValidator.htm
			${WIN_PKG_ROOT}/npp/LINEDRAW.TTF
		)
		IF(MINGW)
			SET(WIN_PKG_FILES ${WIN_PKG_FILES} ${WIN_PKG_ROOT}/libgmp-3.dll ${WIN_PKG_ROOT}/libgmpxx-4.dll)
		ELSE(MINGW)
			SET(WIN_PKG_FILES ${WIN_PKG_FILES} ${WIN_PKG_ROOT}/gmp.dll)
		ENDIF(MINGW)
		INSTALL(FILES ${WIN_PKG_FILES} DESTINATION .)
		IF(NOT MINGW)
			SET(WIN_MSVC_DLL_PATH ${CMAKE_SOURCE_DIR}/windows_package/Microsoft.VC90.CRT)
			SET(WIN_MSVC_DLL_FILES
				${WIN_MSVC_DLL_PATH}/Microsoft.VC90.CRT.manifest
				${WIN_MSVC_DLL_PATH}/msvcm90.dll
				${WIN_MSVC_DLL_PATH}/msvcp90.dll
				${WIN_MSVC_DLL_PATH}/msvcr90.dll
			)
			INSTALL(FILES ${WIN_MSVC_DLL_FILES} DESTINATION Microsoft.VC90.CRT)
		ENDIF(NOT MINGW)
	ENDIF(WIN32)
	# Mandatory linking.
	SET(MANDATORY_LIBRARIES ${MANDATORY_LIBRARIES} piranha)
	ADD_SUBDIRECTORY(pyranha)
ENDIF(BUILD_PYRANHA)

# Install phase: install theories of motion.
SET(THEORIES_FILES
	theories_of_motion/elp1.dps
	theories_of_motion/elp2.dfs
	theories_of_motion/elp3.dfs
)
INSTALL(FILES ${THEORIES_FILES} DESTINATION ${THEORIES_INSTALL_PATH})

INCLUDE(CPack)
