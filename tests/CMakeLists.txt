SET(CMAKE_SKIP_BUILD_RPATH TRUE)
SET(PIRANHA_TESTS_DIRECTORY ${CMAKE_HOME_DIRECTORY}/tests)
#
#
#
#SET(CMAKE_CTEST_COMMAND ctest -V)
#ADD_CUSTOM_TARGET(RUN_TEST_PREPARE COMMAND ${CMAKE_CTEST_COMMAND})
#
# automatic test ecution using ctest unde VS doesn't seem to work properly . Working deirctory!
#
ADD_EXECUTABLE(test_fateman_1 fateman_poly_test1.cpp)
TARGET_INCLUDE_DIRECTORIES(test_fateman_1 PRIVATE ${CMAKE_HOME_DIRECTORY}/src)
TARGET_LINK_LIBRARIES(test_fateman_1 ${MANDATORY_LIBRARIES} piranha_static)

ADD_EXECUTABLE(test_fateman_2 fateman_poly_test2.cpp)
TARGET_INCLUDE_DIRECTORIES(test_fateman_2 PRIVATE ${CMAKE_HOME_DIRECTORY}/src)
TARGET_LINK_LIBRARIES(test_fateman_2 ${MANDATORY_LIBRARIES} piranha_static)

#ADD_EXECUTABLE(small_poly_test small_poly_test.cpp)
#TARGET_LINK_LIBRARIES(small_poly_test ${MANDATORY_LIBRARIES} piranha_static)

ADD_EXECUTABLE(test_gastineau gastineau_fs_test1.cpp)
TARGET_INCLUDE_DIRECTORIES(test_gastineau PRIVATE ${CMAKE_HOME_DIRECTORY}/src)
TARGET_LINK_LIBRARIES(test_gastineau ${MANDATORY_LIBRARIES} piranha_static)

ADD_EXECUTABLE(test_pearce_1  pearce_sparse_poly_test1.cpp)
TARGET_INCLUDE_DIRECTORIES(test_pearce_1 PRIVATE ${CMAKE_HOME_DIRECTORY}/src)
TARGET_LINK_LIBRARIES(test_pearce_1 ${MANDATORY_LIBRARIES} piranha_static)

ADD_EXECUTABLE(test_pearce_2 pearce_sparse_poly_test2.cpp)
TARGET_INCLUDE_DIRECTORIES(test_pearce_2 PRIVATE ${CMAKE_HOME_DIRECTORY}/src)
TARGET_LINK_LIBRARIES(test_pearce_2 ${MANDATORY_LIBRARIES} piranha_static)

ADD_EXECUTABLE(test_celmec celmec_test_01.cpp)
TARGET_INCLUDE_DIRECTORIES(test_celmec PRIVATE ${CMAKE_HOME_DIRECTORY}/src)
TARGET_LINK_LIBRARIES(test_celmec ${MANDATORY_LIBRARIES} piranha_static)

ADD_EXECUTABLE(test_binomial binomial_20000.cpp)
TARGET_INCLUDE_DIRECTORIES(test_binomial PRIVATE ${CMAKE_HOME_DIRECTORY}/src)
TARGET_LINK_LIBRARIES(test_binomial ${MANDATORY_LIBRARIES} piranha_static)

ADD_EXECUTABLE(test_sub EXCLUDE_FROM_ALL sub_tests.cpp)
TARGET_INCLUDE_DIRECTORIES(test_sub PRIVATE ${CMAKE_HOME_DIRECTORY}/src)
TARGET_LINK_LIBRARIES(test_sub ${MANDATORY_LIBRARIES} piranha_static)

ADD_EXECUTABLE(test_sf_01 sf_tests_01.cpp)
TARGET_INCLUDE_DIRECTORIES(test_sf_01 PRIVATE ${CMAKE_HOME_DIRECTORY}/src)
TARGET_LINK_LIBRARIES(test_sf_01 ${MANDATORY_LIBRARIES} piranha_static)

ADD_EXECUTABLE(test_dfs1 dfs_mult_01.cpp)
TARGET_INCLUDE_DIRECTORIES(test_dfs1 PRIVATE ${CMAKE_HOME_DIRECTORY}/src)
TARGET_LINK_LIBRARIES(test_dfs1 ${MANDATORY_LIBRARIES} piranha_static)

ADD_EXECUTABLE(test_dfs2 dfs_mult_02.cpp)
TARGET_INCLUDE_DIRECTORIES(test_dfs2 PRIVATE ${CMAKE_HOME_DIRECTORY}/src)
TARGET_LINK_LIBRARIES(test_dfs2 ${MANDATORY_LIBRARIES} piranha_static)

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

#set_property(TEST test_gastineau PROPERTY WORKING_DIRECTORY ${PIRANHA_TESTS_DIRECTORY})
