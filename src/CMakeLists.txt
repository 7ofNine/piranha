#
#src/CMakeList
#
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/template.settings.cpp ${CMAKE_CURRENT_SOURCE_DIR}/settings.cpp) # generate the settings.cpp file from the detected parameters

SET(SRC_LIST
        degree_truncator.cpp
        norm_truncator.cpp
        math.cpp
        stats.cpp
        settings.cpp
    )

SET(SRC_LIST_HEADERS
        core/mp.h core/psym.h core/truncators/degree.h
        core/exceptions.h core/settings.h core/truncators/norm.h
        core/exceptions.h core/mp.h
        core/base_classes/base_counting_allocator.h
        core/config.h core/exceptions.h core/memory.h core/settings.h
        core/Psym.h
        core/stats.h
)

SET(SRC_TEMPLATE template.settings.cpp)
SET_SOURCE_FILES_PROPERTIES(template.settings.cpp PROPERTIES HEADER_FILE_ONLY ON) #disable compilation. Only used for generation
SET(SRC_LIST ${SRC_LIST} ${SRC_LIST_HEADERS} ${SRC_TEMPLATE})  # combine src and headers and template(s) for MSVC projects

SOURCE_GROUP(TREE ${CMAKE_CURRENT_LIST_DIR} PREFIX "Templates" FILES ${SRC_TEMPLATE})



# define target piranha.dll
#
ADD_COMPILE_OPTIONS(/wd4251)   # suppress dll interface warnings for stl items. 

ADD_LIBRARY(piranha SHARED ${SRC_LIST} ${SRC_LIST_HEADERS}) # define target piranha.dll
set_property(TARGET piranha PROPERTY CXX_STANDARD_REQUIRED YES)

#TARGET_INCLUDE_DIRECTORIES(piranha PRIVATE ${Boost_INCLUDE_DIR})  # see below TARGET_LINK_LIBRARIES
SET_PROPERTY(TARGET piranha PROPERTY COMPILE_DEFINITIONS PIRANHA_DLL_EXPORT_API) # The dll exports its intefaces
TARGET_LINK_LIBRARIES(piranha ${MANDATORY_LIBRARIES} Boost::thread Boost::disable_autolinking Boost::dynamic_linking)  #this seems also to set the boost include directory. Is that specified somewhere ?



# copy libraries so that we actually can run the tests and/or the single projects
if(BUILD_TESTS)
    # somehow ENVIRONMENT for tests properties doesn't wotk (yet?)
    # CMAKE still doesn't support copying runtime dependencies properly to the targets. It's a pain
    get_filename_component(BASE_NAME_THREAD ${Boost_THREAD_LIBRARY_RELEASE} NAME_WE)

    add_custom_command(TARGET piranha POST_BUILD
                        COMMAND ${CMAKE_COMMAND} -E copy_if_different 
                                "$<TARGET_FILE:piranha>"
                                "${CMAKE_BINARY_DIR}/tests/$<CONFIG>/$<TARGET_FILE_NAME:piranha>"
                        COMMAND ${CMAKE_COMMAND} -E copy_if_different
                                "${MPIR_INCLUDE_DIR}/mpir.dll" "${CMAKE_BINARY_DIR}/tests/$<CONFIG>/mpir.dll" 
                        COMMAND ${CMAKE_COMMAND} -E copy_if_different 
                                "${Boost_LIBRARY_DIR_RELEASE}/${BASE_NAME_THREAD}.dll"
                                "${CMAKE_BINARY_DIR}/tests/$<CONFIG>/${BASE_NAME_THREAD}.dll"
    ) 
endif()


if(BUILD_PYRANHA)
    get_filename_component(BASE_NAME_THREAD ${Boost_THREAD_LIBRARY_RELEASE} NAME_WE)
    
    string(TOUPPER ${PIRANHA_PYTHON_BASE} PIRANHA_PYTHON_BASE_UPPER)
    get_filename_component(BASE_NAME_PYTHON ${Boost_${PIRANHA_PYTHON_BASE_UPPER}_LIBRARY_RELEASE} NAME_WE)

    INSTALL(TARGETS piranha RUNTIME DESTINATION .)
    SET(PYRANA_BOOST_INSTALL_FILES "${Boost_LIBRARY_DIR_RELEASE}/${BASE_NAME_THREAD}.dll" "${Boost_LIBRARY_DIR_RELEASE}/${BASE_NAME_PYTHON}.dll" )
    INSTALL(FILES  ${PYRANA_BOOST_INSTALL_FILES} DESTINATION . )
endif()
      

# Target to get the header files into MSVC      
SET(INCS
        piranha.h
        core/utils.h
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
        core/lambdas.h
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

SET(PYRANHA_PYTHON
    ../pyranha/detail.py
    ../pyranha/global.__init__.py
    ../pyranha/template.__init__.py
    ../pyranha/template.global.__init__.py
    ../pyranha/Celmec/__init__.py
    ../pyranha/Celmec/detail.py
    ../pyranha/Celmec/low_thrust.py
    ../pyranha/Test/__init__.py
    ../pyranha/Test/detail.py
    ../extras/template.ipy_profile_pyranha.py
) 
SOURCE_GROUP(TREE ${CMAKE_SOURCE_DIR} PREFIX "Python" FILES ${PYRANHA_PYTHON})   

ADD_LIBRARY(Headers STATIC ${INCS} ${PYRANHA_PYTHON})
TARGET_INCLUDE_DIRECTORIES(Headers PRIVATE ${Boost_INCLUDE_DIR})
SET_TARGET_PROPERTIES(Headers PROPERTIES LINKER_LANGUAGE CXX)
