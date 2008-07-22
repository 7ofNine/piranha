SET(SRC_LIST "")

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/template.settings.cpp ${CMAKE_CURRENT_SOURCE_DIR}/settings.cpp)

SET(SRC_LIST ${SRC_LIST}
	shared_args.cpp
	base_degree_truncator.cpp
	base_expo_truncator.cpp
	base_norm_truncator.cpp
	settings.cpp
	stream_manager.cpp
	psym.cpp
)

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
		tbb/recursive_mutex.cpp
		tbb/spin_mutex.cpp
		tbb/spin_rw_mutex.cpp
		tbb/task.cpp
		tbb/tbb_misc.cpp
		tbb/tbb_thread.cpp
	)
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
		boost_python/object/function_doc_signature.cpp
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
	ADD_DEFINITIONS(-DBOOST_PYTHON_SOURCE)
	ADD_LIBRARY(piranha SHARED ${SRC_LIST} ${SHARED_SRC_LIST})
	SET_TARGET_PROPERTIES(piranha PROPERTIES LINKER_LANGUAGE CXX)
	TARGET_LINK_LIBRARIES(piranha ${MANDATORY_LIBRARIES} ${PYTHON_LIBRARY})
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
