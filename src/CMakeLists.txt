SET(SRC_LIST "")
IF(NOT GCC_CXX_COMPILER AND NOT ICC_CXX_COMPILER AND BUILD_LNP)
  MESSAGE(FATAL_ERROR "The C++ STL extension 'ext/slist' is needed to compile lnp series.")
ELSE(NOT GCC_CXX_COMPILER AND NOT ICC_CXX_COMPILER AND BUILD_LNP)
  MESSAGE(STATUS "Header 'ext/slist' is available.")
ENDIF(NOT GCC_CXX_COMPILER AND NOT ICC_CXX_COMPILER AND BUILD_LNP)
IF(BUILD_TASS)
  IF(BUILD_FORTRAN)
    CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/tass17/posired7.f.template
      ${CMAKE_CURRENT_SOURCE_DIR}/tass17/posired7.f)
    SET(SRC_LIST ${SRC_LIST} tass17/posired7.f)
  ENDIF(BUILD_FORTRAN)
  SET(SRC_LIST ${SRC_LIST} tass17/tass17.cpp)
ENDIF(BUILD_TASS)
IF(BUILD_KEPLERIAN_PROCESSOR)
  SET(SRC_LIST ${SRC_LIST} keplerian_processor/keplerian_processor.cpp)
ENDIF(BUILD_KEPLERIAN_PROCESSOR)
IF(BUILD_SP)
  SET(SRC_LIST ${SRC_LIST} symbol_limiter.cpp)
ENDIF(BUILD_SP)
# Order is important here: in tass.cpp, for instance, we have destruction of series that depends on
# psymbol_manager to be present. Hence psymbol.cpp must go after...
SET(SRC_LIST ${SRC_LIST} settings_manager.cpp stats.cpp astro.cpp arg_manager.cpp
  stream_manager.cpp psymbol.cpp)

IF(BUILD_TBB_MULTITHREADING)
  SET(SRC_LIST ${SRC_LIST}
    tbb/cache_aligned_allocator.cpp
    tbb/concurrent_hash_map.cpp
    tbb/concurrent_queue.cpp
    tbb/concurrent_vector.cpp
    tbb/itt_notify.cpp
    tbb/mutex.cpp
    tbb/pipeline.cpp
    tbb/queuing_mutex.cpp
    tbb/queuing_rw_mutex.cpp
    tbb/spin_mutex.cpp
    tbb/spin_rw_mutex.cpp
    tbb/task.cpp
    tbb/tbb_misc.cpp
  )
  IF(WIN32)
    ADD_DEFINITIONS(-DUSE_WINTHREAD)
  ELSE(WIN32)
    ADD_DEFINITIONS(-DUSE_PTHREAD)
    SET(MANDATORY_LIBRARIES ${MANDATORY_LIBRARIES} pthread dl)
  ENDIF(WIN32)
ENDIF(BUILD_TBB_MULTITHREADING)

IF(BUILD_PYRANHA)
  SET(SHARED_SRC_LIST
    boost_python/converter/arg_to_python_base.cpp
    boost_python/converter/builtin_converters.cpp
    boost_python/converter/from_python.cpp
    boost_python/converter/registry.cpp
    boost_python/converter/type_id.cpp
    boost_python/object/class.cpp
    boost_python/object/enum.cpp
    boost_python/object/function.cpp
    boost_python/object/inheritance.cpp
    boost_python/object/iterator.cpp
    boost_python/object/life_support.cpp
    boost_python/object/pickle_support.cpp
    boost_python/object/stl_iterator.cpp
    boost_python/dict.cpp
    boost_python/errors.cpp
    boost_python/exec.cpp
    boost_python/import.cpp
    boost_python/list.cpp
    boost_python/long.cpp
    boost_python/module.cpp
    boost_python/numeric.cpp
    boost_python/object_operators.cpp
    boost_python/object_protocol.cpp
    boost_python/slice.cpp
    boost_python/str.cpp
    boost_python/tuple.cpp
    boost_python/wrapper.cpp
    )
  ADD_LIBRARY(piranha SHARED ${SRC_LIST} ${SHARED_SRC_LIST})
  SET_TARGET_PROPERTIES(piranha PROPERTIES LINKER_LANGUAGE CXX)
  TARGET_LINK_LIBRARIES(piranha ${MANDATORY_LIBRARIES})
  INSTALL(TARGETS piranha
    RUNTIME DESTINATION ${LIB_INSTALL_PATH}
    LIBRARY DESTINATION ${LIB_INSTALL_PATH}
  )
ENDIF(BUILD_PYRANHA)

IF(BUILD_MAIN)
  ADD_LIBRARY(piranha_static STATIC ${SRC_LIST})
  TARGET_LINK_LIBRARIES(piranha_static ${MANDATORY_LIBRARIES})
  SET_TARGET_PROPERTIES(piranha_static PROPERTIES LINKER_LANGUAGE CXX)
ENDIF(BUILD_MAIN)
