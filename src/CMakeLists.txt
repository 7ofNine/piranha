SET(SRC_LIST "")

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/template.settings.cpp ${CMAKE_CURRENT_SOURCE_DIR}/settings.cpp)

SET(SRC_LIST ${SRC_LIST}
	degree_truncator.cpp
	norm_truncator.cpp
	math.cpp
	memory.cpp
	settings.cpp
	Psym.cpp
	stats.cpp
)

#ADD_DEFINITIONS(-DBOOST_THREAD_BUILD_LIB)
IF(WIN32)
	SET(SRC_LIST ${SRC_LIST}
		boost_thread/win32/thread.cpp
		boost_thread/win32/tss_dll.cpp
		boost_thread/win32/tss_pe.cpp
	)
ELSE(WIN32)
	SET(SRC_LIST ${SRC_LIST}
		boost_thread/pthread/once.cpp
		boost_thread/pthread/thread.cpp
	)
ENDIF(WIN32)

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
	ADD_DEFINITIONS(-D_PIRANHA_DLL_EXPORT_API -DBOOST_PYTHON_SOURCE)
	ADD_LIBRARY(piranha SHARED ${SRC_LIST} ${SHARED_SRC_LIST})
	SET_PROPERTY(TARGET piranha PROPERTY COMPILE_DEFINITIONS BOOST_THREAD_BUILD_DLL)
	SET_TARGET_PROPERTIES(piranha PROPERTIES LINKER_LANGUAGE CXX)
	TARGET_LINK_LIBRARIES(piranha ${MANDATORY_LIBRARIES} ${PYTHON_LIBRARY})
	INSTALL(TARGETS piranha
		RUNTIME DESTINATION ${LIB_INSTALL_PATH}
		LIBRARY DESTINATION ${LIB_INSTALL_PATH}
	)
ENDIF(BUILD_PYRANHA)

IF(BUILD_MAIN OR BUILD_TESTS)
	ADD_LIBRARY(piranha_static STATIC ${SRC_LIST})
	SET_PROPERTY(TARGET piranha_static PROPERTY COMPILE_DEFINITIONS BOOST_THREAD_BUILD_LIB)
	TARGET_LINK_LIBRARIES(piranha_static ${MANDATORY_LIBRARIES})
	SET_TARGET_PROPERTIES(piranha_static PROPERTIES LINKER_LANGUAGE CXX)
ENDIF(BUILD_MAIN OR BUILD_TESTS)

#FILE(GLOB_RECURSE INCS "*.h")

#ADD_CUSTOM_TARGET(Headers SOURCES ${INCS})
SET(INCS
				piranha.h
				core/utils.h
				core/type_traits.h
				core/stats.h
				core/settings.h
				core/Psym.h
				core/power_cache.h
				core/p_exceptions.h
				core/ntuple.h
				core/mp.h
				core/memory.h
				core/math.h
				core/integer_typedefs.h
				core/exceptions.h
				core/config.h
				core/common_functors.h
				core/coded_hash_table.h
				core/truncators/degree.h
				core/truncators/norm.h
				core/truncators/power_series.h
				core/truncators/truncators.h
				core/polynomial/polynomial.h
				core/polynomial/named_polynomial.h
				core/polynomial_common/base_polynomial.h
				core/polynomial_common/expo_vector.h
				core/polynomial_common/expo_vector_mp.h
				core/polynomial_common/monomial.h
				core/polynomial_common/polynomial_multiplier.h
				core/polynomial_cf/common_polynomial_cf_toolbox.h
				core/polynomial_cf/polynomial_cf.h
				core/poisson_series_common/celmec_toolbox.h
				core/poisson_series_common/common_poisson_series_toolbox.h
				core/poisson_series_common/jacobi_anger_toolbox.h
				core/poisson_series_common/poisson_series_multiplier.h
				core/poisson_series/poisson_series.h
				core/numerical_coefficients/double_cf.h
				core/numerical_coefficients/mpf_cf.h
				core/numerical_coefficients/mpq_cf.h
				core/numerical_coefficients/mpz_cf.h
				core/harmonic_series/base_harmonic_series.h
				core/harmonic_series/named_harmonic_series.h
				core/harmonic_series/trig_vector.h
				core/harmonic_series/trig_vector_mp.h
				core/fourier_series/base_fourier_series.h
				core/fourier_series/common_fourier_series_toolbox.h
				core/fourier_series/fourier_series.h
				core/fourier_series/fourier_series_term.h
				core/fourier_series/named_fourier_series.h
				core/base_classes/base_series_complex_toolbox.h
				core/base_classes/base_series_def.h
				core/base_classes/base_series_io.h
				core/base_classes/base_series_manip.h
				core/base_classes/base_series_math.h
				core/base_classes/base_series_mp.h
				core/base_classes/base_series_multiplier.h
				core/base_classes/base_series_multiplier_mp.h
				core/base_classes/base_series_probe.h
				core/base_classes/base_series_special_functions.h
				core/base_classes/base_series_tag.h
				core/base_classes/base_term.h
				core/base_classes/binomial_exponentiation_toolbox.h
				core/base_classes/cf_power_series.h
				core/base_classes/cf_series.h
				core/base_classes/cf_series_complex_toolbox.h
				core/base_classes/cf_series_io.h
				core/base_classes/cf_series_manip.h
				core/base_classes/cf_series_math.h
				core/base_classes/cf_series_probe.h
				core/base_classes/cf_series_special_functions.h
				core/base_classes/coded_multiplier.h
				core/base_classes/coded_multiplier_mp.h
				core/base_classes/common_args_descriptions.h
				core/base_classes/named_power_series.h
				core/base_classes/named_series.h
				core/base_classes/named_series_complex_toolbox.h
				core/base_classes/named_series_def.h
				core/base_classes/named_series_io.h
				core/base_classes/named_series_manip.h
				core/base_classes/named_series_math.h
				core/base_classes/named_series_mp.h
				core/base_classes/named_series_probe.h
				core/base_classes/named_series_special_functions.h
				core/base_classes/null_truncator.h
				core/base_classes/numerical_container.h
				core/base_classes/numerical_container_complex_toolbox.h
				core/base_classes/numerical_container_mp.h
				core/base_classes/numerical_container_tag.h
				core/base_classes/series_multiplication.h
				core/base_classes/series_multiplier.h
				core/base_classes/vector_key.h
				core/atomic_counters/atomic_counter_gcc_41.h
				core/atomic_counters/atomic_counter_generic.h
				core/atomic_counters/atomic_counter_msvc_long.h
				core/atomic_counters/atomic_counters.h
				core/mp/complex_generic_mp_container.h
				core/mp/mp_commons.h
				core/mp/piranha_gmp.h
				manipulators/dfs.h
				manipulators/dpoly.h
				manipulators/dps.h
				manipulators/dqps.h
				manipulators/qpoly.h
				manipulators/qps.h
				manipulators/qqpoly.h
				manipulators/qqps.h
				manipulators/zpoly.h
				../pyranha/args_tuple.h
				../pyranha/boost_python_container_conversions.h
				../pyranha/boost_python_p_exceptions.h
				../pyranha/commons.h
				../pyranha/exceptions.h
				../pyranha/mp_classes.h
				../pyranha/series_instantiations.h
				../pyranha/stl_containers.h
				)
				
ADD_CUSTOM_TARGET(Headers SOURCES ${INCS})