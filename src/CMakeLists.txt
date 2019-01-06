SET(SRC_LIST "")

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/settings.cpp.template ${CMAKE_CURRENT_SOURCE_DIR}/settings.cpp) # generate the settings.cpp file from the detected paramters

SET(SRC_LIST ${SRC_LIST} degree_truncator.cpp
                         norm_truncator.cpp
                         math.cpp
                         memory.cpp
                         settings.cpp
                         Psym.cpp
                         stats.cpp
                         settings.cpp.template
)

SET(SRC_LIST_HEADERS core/mp.h core/psym.h core/truncators/degree.h
                     core/exceptions.h core/settings.h core/truncators/norm.h
                     core/exceptions.h core/mp.h
                     core/atomic_counters/atomic_counters.h core/base_classes/base_counting_allocator.h
                     core/config.h core/exceptions.h core/memory.h core/settings.h
                     core/Psym.h
                     core/stats.h
)

SET(SRC_LIST ${SRC_LIST} ${SRC_LIST_HEADERS})  # combine src and headers for MSVC projects

# find Boost
SET(Boost_USE_RELEASE_LIBS ON)
SET(Boost_USE_STATIC_LIBS OFF)
SET(Boost_USE_STATIC_RUNTIME OFF)
set(Boost_USE_MULTITHREADED ON)
FIND_PACKAGE(Boost 1.66.0 REQUIRED COMPONENTS thread python27)
if(Boost_FOUND)
	MESSAGE(STATUS "Found Boost libraries: ${Boost_LIBRARIES}")
else()
	MESSAGE(FATAL "Couldn't find Boost libraries")
endif()


# define target piranha.dll
ADD_LIBRARY(piranha SHARED ${SRC_LIST} ${SRC_LIST_HEADERS}) # define target piranha.dll

TARGET_INCLUDE_DIRECTORIES(piranha PRIVATE ${PYTHON_INCLUDE_DIR} ${Boost_INCLUDE_DIR})
SET_PROPERTY(TARGET piranha PROPERTY COMPILE_DEFINITIONS _PIRANHA_DLL_EXPORT_API) # The dll exports its intefraces
SET_TARGET_PROPERTIES(piranha PROPERTIES LINKER_LANGUAGE CXX CXX_STANDARD 11)
TARGET_LINK_LIBRARIES(piranha ${MANDATORY_LIBRARIES} ${PYTHON_LIBRARY} Boost::thread Boost::python27)

# copy libraries so that we actually can run the tests and/or the single projects
if(BUILD_TESTS)
# CMAKE still doesn't support copying runtime dependencies properly to the targets. It's a pain
get_filename_component(BASE_NAME_THREAD ${Boost_THREAD_LIBRARY_RELEASE} NAME_WE)
get_filename_component(BASE_NAME_PYTHON ${Boost_PYTHON27_LIBRARY_RELEASE} NAME_WE)
#get_filename_component(BASE_NAME_SYSTEM ${Boost_SYSTEM_LIBRARY_RELEASE} NAME_WE)

	add_custom_command(TARGET piranha POST_BUILD
						COMMAND ${CMAKE_COMMAND} -E copy_if_different 
						        "$<TARGET_FILE:piranha>"
								"${CMAKE_BINARY_DIR}/tests/$<CONFIG>/$<TARGET_FILE_NAME:piranha>"
						COMMAND ${CMAKE_COMMAND} -E copy_if_different
								"${MPIR_INCLUDE_DIR}/mpir.dll" "${CMAKE_BINARY_DIR}/tests/$<CONFIG>/mpir.dll" 
						COMMAND ${CMAKE_COMMAND} -E copy_if_different 
								"${Boost_LIBRARY_DIR_RELEASE}/${BASE_NAME_THREAD}.dll"
								"${CMAKE_BINARY_DIR}/tests/$<CONFIG>/${BASE_NAME_THREAD}.dll"
						COMMAND ${CMAKE_COMMAND} -E copy_if_different 
								"${Boost_LIBRARY_DIR_RELEASE}/${BASE_NAME_PYTHON}.dll"
								"${CMAKE_BINARY_DIR}/tests/$<CONFIG>/${BASE_NAME_PYTHON}.dll"
								# boost:system is no longer a dll in Boost >= 1.69
#						COMMAND ${CMAKE_COMMAND} -E copy_if_different 
#								"${Boost_LIBRARY_DIR_RELEASE}/${BASE_NAME_SYSTEM}.dll"
#								"${CMAKE_BINARY_DIR}/tests/$<CONFIG>/${BASE_NAME_SYSTEM}.dll"			
#								
					  )	
endif()

	
INSTALL(TARGETS piranha 
				RUNTIME DESTINATION . 
                LIBRARY DESTINATION . 
		)

# Target to get the header files into MSVC		
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
                core/base_classes/base_counting_allocator.h
                core/base_classes/base_power_series.h
                core/base_classes/base_series.h
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
SOURCE_GROUP(TREE ${CMAKE_SOURCE_DIR} PREFIX "Header Files" FILES ${INCS})
ADD_LIBRARY(Headers STATIC ${INCS})
TARGET_INCLUDE_DIRECTORIES(Headers PRIVATE ${Boost_INCLUDE_DIR})
SET_TARGET_PROPERTIES(Headers PROPERTIES LINKER_LANGUAGE CXX CXX_STANDARD 11)
