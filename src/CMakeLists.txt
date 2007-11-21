SET(SRC_LIST "")
IF(BUILD_FORTRAN)
  CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/tass17/posired7.f.template
    ${CMAKE_CURRENT_SOURCE_DIR}/tass17/posired7.f)
  SET(SRC_LIST ${SRC_LIST} tass17/posired7.f)
ENDIF(BUILD_FORTRAN)
IF(BUILD_PYRANHA_SP)
  SET(SRC_LIST ${SRC_LIST} symbol_limiter.cpp)
ENDIF(BUILD_PYRANHA_SP)
# Order is important here: in tass.cpp, for instance, we have destruction of series that depends on
# psymbol_manager to be present. Hence psymbol.cpp must go after...
SET(SRC_LIST ${SRC_LIST}
  norm_truncation.cpp
  buffer.cpp
  settings_manager.cpp
  stats.cpp
  astro.cpp
  arg_manager.cpp
  stream_manager.cpp
  psymbol.cpp
)

IF(ENABLE_HOARD)
  INCLUDE_DIRECTORIES("hoard")
  INCLUDE_DIRECTORIES("hoard/heaplayers")
  INCLUDE_DIRECTORIES("hoard/heaplayers/util")
  SET(MANDATORY_LIBRARIES ${MANDATORY_LIBRARIES} pthread dl)
  SET(SRC_LIST ${SRC_LIST} hoard/libhoard.cpp)
ENDIF(ENABLE_HOARD)

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
    piranha_tbb.cpp
  )
  IF(WIN32)
# This is irrelevant at the moment, since TBB does not compile yet in MinGW.
    ADD_DEFINITIONS(-DUSE_WINTHREAD)
  ELSE(WIN32)
    ADD_DEFINITIONS(-DUSE_PTHREAD)
    SET(MANDATORY_LIBRARIES ${MANDATORY_LIBRARIES} pthread dl)
  ENDIF(WIN32)
  ADD_DEFINITIONS(-DPIRANHA_TBB)
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
  #ADD_DEFINITIONS(-DBOOST_PYTHON_DYNAMIC_LIB)
  ADD_LIBRARY(piranha SHARED ${SRC_LIST} ${SHARED_SRC_LIST})
  SET_TARGET_PROPERTIES(piranha PROPERTIES LINKER_LANGUAGE CXX)
  TARGET_LINK_LIBRARIES(piranha ${MANDATORY_LIBRARIES})
  INSTALL(TARGETS piranha
    RUNTIME DESTINATION ${LIB_INSTALL_PATH}
    LIBRARY DESTINATION ${LIB_INSTALL_PATH}
  )
ENDIF(BUILD_PYRANHA)

IF(BUILD_MAIN OR BUILD_TESTS)
  ADD_LIBRARY(piranha_static STATIC ${SRC_LIST})
  TARGET_LINK_LIBRARIES(piranha_static ${MANDATORY_LIBRARIES})
  SET_TARGET_PROPERTIES(piranha_static PROPERTIES LINKER_LANGUAGE CXX)
ENDIF(BUILD_MAIN OR BUILD_TESTS)
