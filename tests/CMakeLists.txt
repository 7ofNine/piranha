SET(CMAKE_SKIP_BUILD_RPATH TRUE)
SET(PIRANHA_TESTS_DIRECTORY ${CMAKE_HOME_DIRECTORY}/tests)
ADD_COMPILE_OPTIONS(/bigobj)
ADD_DEFINITIONS(-D_PIRANHA_DLL_IMPORT_API) # we import piranha.dll
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
list(APPEND testincludePath "${CMAKE_HOME_DIRECTORY}/src" "${CMAKE_HOME_DIRECTORY}/src/core")
#
#
#
#SET(CMAKE_CTEST_COMMAND ctest -V)
#ADD_CUSTOM_TARGET(RUN_TEST_PREPARE COMMAND ${CMAKE_CTEST_COMMAND})
#
# automatic test execution using ctest unde VS doesn't seem to work properly
# we have to copy somehow the mpir.dll and mpfr.dll into the DEBUG/RELEASE directories for the tests to run. Automate!
#
# Temporarily turn off the old tests until complete tests are possible
# they are nice examples but no complete test

#ADD_EXECUTABLE(test_fateman_1 fateman_poly_test1.cpp)
#TARGET_INCLUDE_DIRECTORIES(test_fateman_1 PRIVATE ${CMAKE_HOME_DIRECTORY}/src)
#TARGET_LINK_LIBRARIES(test_fateman_1 piranha)

#ADD_EXECUTABLE(test_fateman_2 fateman_poly_test2.cpp)
#TARGET_INCLUDE_DIRECTORIES(test_fateman_2 PRIVATE ${CMAKE_HOME_DIRECTORY}/src)
#TARGET_LINK_LIBRARIES(test_fateman_2 piranha)

#ADD_EXECUTABLE(small_poly_test small_poly_test.cpp)
#TARGET_LINK_LIBRARIES(small_poly_test ${MANDATORY_LIBRARIES} piranha)

#ADD_EXECUTABLE(test_gastineau gastineau_fs_test1.cpp)
#TARGET_INCLUDE_DIRECTORIES(test_gastineau PRIVATE ${CMAKE_HOME_DIRECTORY}/src)
#TARGET_LINK_LIBRARIES(test_gastineau piranha)

#ADD_EXECUTABLE(test_pearce_1  pearce_sparse_poly_test1.cpp)
#TARGET_INCLUDE_DIRECTORIES(test_pearce_1 PRIVATE ${CMAKE_HOME_DIRECTORY}/src)
#TARGET_LINK_LIBRARIES(test_pearce_1 piranha)

#ADD_EXECUTABLE(test_pearce_2 pearce_sparse_poly_test2.cpp)
#TARGET_INCLUDE_DIRECTORIES(test_pearce_2 PRIVATE ${CMAKE_HOME_DIRECTORY}/src)
#TARGET_LINK_LIBRARIES(test_pearce_2 piranha)

ADD_EXECUTABLE(test_celmec celmec_test_01.cpp)
TARGET_INCLUDE_DIRECTORIES(test_celmec PRIVATE ${CMAKE_HOME_DIRECTORY}/src)
TARGET_LINK_LIBRARIES(test_celmec piranha)

#ADD_EXECUTABLE(test_binomial binomial_20000.cpp)
#TARGET_INCLUDE_DIRECTORIES(test_binomial PRIVATE ${CMAKE_HOME_DIRECTORY}/src)
#TARGET_LINK_LIBRARIES(test_binomial piranha)

ADD_EXECUTABLE(test_sub sub_tests.cpp)
TARGET_INCLUDE_DIRECTORIES(test_sub PRIVATE ${CMAKE_HOME_DIRECTORY}/src)
TARGET_LINK_LIBRARIES(test_sub piranha)

#ADD_EXECUTABLE(test_sf_01 sf_tests_01.cpp)
#TARGET_INCLUDE_DIRECTORIES(test_sf_01 PRIVATE ${CMAKE_HOME_DIRECTORY}/src)
#TARGET_LINK_LIBRARIES(test_sf_01 piranha)

#ADD_EXECUTABLE(test_dfs1 dfs_mult_01.cpp)
#TARGET_INCLUDE_DIRECTORIES(test_dfs1 PRIVATE ${CMAKE_HOME_DIRECTORY}/src)
#TARGET_LINK_LIBRARIES(test_dfs1 piranha)

#ADD_EXECUTABLE(test_dfs2 dfs_mult_02.cpp)
#TARGET_INCLUDE_DIRECTORIES(test_dfs2 PRIVATE ${CMAKE_HOME_DIRECTORY}/src)
#TARGET_LINK_LIBRARIES(test_dfs2 piranha)

ADD_EXECUTABLE(test_read test_read.cpp)
TARGET_INCLUDE_DIRECTORIES(test_read PRIVATE ${CMAKE_HOME_DIRECTORY}/src)
TARGET_LINK_LIBRARIES(test_read piranha)


MACRO(CREATE_BOOST_TEST arg1)
	ADD_EXECUTABLE(${arg1} ${arg1}.cpp)
	TARGET_INCLUDE_DIRECTORIES(${arg1} PRIVATE ${CMAKE_HOME_DIRECTORY}/src)
	TARGET_LINK_LIBRARIES(${arg1} piranha)
	ADD_TEST(${arg1} ${arg1})
	SET_TARGET_PROPERTIES(${arg1} PROPERTIES FOLDER Tests)
ENDMACRO(CREATE_BOOST_TEST)


## new boost tests
#ADD_EXECUTABLE(test_psym test_psym.cpp)
#TARGET_INCLUDE_DIRECTORIES(test_psym PRIVATE ${CMAKE_HOME_DIRECTORY}/src)
#TARGET_LINK_LIBRARIES(test_psym piranha)
#ADD_TEST(test_psym test_psym)

#ADD_EXECUTABLE(test_vector_key test_vector_key.cpp)
#TARGET_INCLUDE_DIRECTORIES(test_vector_key PRIVATE ${CMAKE_HOME_DIRECTORY}/src)
#TARGET_LINK_LIBRARIES(test_vector_key piranha)
#ADD_TEST(test_vector_key test_vector_key)

#ADD_EXECUTABLE(test_expo_vector test_expo_vector.cpp)
#TARGET_INCLUDE_DIRECTORIES(test_expo_vector PRIVATE ${CMAKE_HOME_DIRECTORY}/src)
#TARGET_LINK_LIBRARIES(test_expo_vector piranha)
#ADD_TEST(test_expo_vector test_expo_vector)

#ADD_EXECUTABLE(test_trig_vector test_trig_vector.cpp)
#TARGET_INCLUDE_DIRECTORIES(test_trig_vector PRIVATE ${CMAKE_HOME_DIRECTORY}/src)
#TARGET_LINK_LIBRARIES(test_trig_vector piranha)
#ADD_TEST(test_trig_vector test_trig_vector)

#ADD_EXECUTABLE(test_base_term test_base_term.cpp)
#TARGET_INCLUDE_DIRECTORIES(test_base_term PRIVATE ${CMAKE_HOME_DIRECTORY}/src)
#TARGET_LINK_LIBRARIES(test_base_term piranha)
#ADD_TEST(test_base_term test_base_term)

CREATE_BOOST_TEST(test_psym)
CREATE_BOOST_TEST(test_vector_key)
CREATE_BOOST_TEST(test_expo_vector)
CREATE_BOOST_TEST(test_trig_vector)
CREATE_BOOST_TEST(test_base_term)
CREATE_BOOST_TEST(test_Monomial)
CREATE_BOOST_TEST(test_BaseSeries)
CREATE_BOOST_TEST(test_Polynomial)
CREATE_BOOST_TEST(test_NamedSeries)
CREATE_BOOST_TEST(test_FourierSeriesTerm)
CREATE_BOOST_TEST(test_FourierSeries)
CREATE_BOOST_TEST(test_PoissonSeries)
CREATE_BOOST_TEST(test_qps_series)
CREATE_BOOST_TEST(test_numericalContainer)
CREATE_BOOST_TEST(test_doubleCf)
# and several more we don't know yet. How to test the specific instantiation qps, dps etc.



# doesn't work . Possibly only for VS express
#set_target_properties( test_fateman_1 PROPERTIES FOLDER TESTS)
# Add tests.
ADD_TEST(test_fateman_1 test_fateman_1)
#ADD_DEPENDENCIES(RUN_TEST_PREPARE test_fateman_1)

ADD_TEST(test_fateman_2 test_fateman_2)

#ADD_TEST(small small_poly_test)  # runs for ever??

#ADD_TEST(test_gastineau COMMAND test_gastineau)

ADD_TEST(test_pearce_1 test_pearce_1 )

#ADD_TEST(test_pearce_2 test_pearce_2)

ADD_TEST(test_binomial test_binomial)
#ADD_DEPENDENCIES(RUN_TEST_PREPARE test_binomial)

ADD_TEST(test_celmec test_celmec)

ADD_TEST(test_sub test_sub )
#ADD_DEPENDENCIES(RUN_TEST_PREPARE test_sub)

ADD_TEST(test_sf_01 test_sf_01)
#ADD_DEPENDENCIES(RUN_TEST_PREPARE test_sf_01)

#ADD_TEST(test_dfs1 test_dfs1 )

#ADD_TEST(test_dfs2 test_dfs2 )

set_tests_properties(test_celmec PROPERTIES ENVIRONMENT "PATH=D:\\Dev\\mpir\\dll\\x64\\Release;$ENV{PATH}") # temp just for testing generation
