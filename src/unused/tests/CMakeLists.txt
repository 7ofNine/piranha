ADD_EXECUTABLE(fateman_test fateman_test.cpp)
TARGET_LINK_LIBRARIES(fateman_test ${MANDATORY_LIBRARIES} piranha_static)

ADD_EXECUTABLE(gastest_poly gastest_poly.cpp)
TARGET_LINK_LIBRARIES(gastest_poly ${MANDATORY_LIBRARIES} piranha_static)

ADD_EXECUTABLE(gastest_elp3 gastest_elp3.cpp)
TARGET_LINK_LIBRARIES(gastest_elp3 ${MANDATORY_LIBRARIES} piranha_static)

ADD_EXECUTABLE(gastest_elp3_dynamic gastest_elp3_dynamic.cpp)
TARGET_LINK_LIBRARIES(gastest_elp3_dynamic ${MANDATORY_LIBRARIES} piranha_static)

ADD_EXECUTABLE(gastest_elp3_ext gastest_elp3_ext.cpp)
TARGET_LINK_LIBRARIES(gastest_elp3_ext ${MANDATORY_LIBRARIES} piranha_static)

ADD_EXECUTABLE(gastest_elp3_dynamic_ext gastest_elp3_dynamic_ext.cpp)
TARGET_LINK_LIBRARIES(gastest_elp3_dynamic_ext ${MANDATORY_LIBRARIES} piranha_static)

ADD_EXECUTABLE(mixed_sin_cos_test mixed_sin_cos_test.cpp)
TARGET_LINK_LIBRARIES(mixed_sin_cos_test ${MANDATORY_LIBRARIES} piranha_static)

ADD_EXECUTABLE(complexp_test complexp_test.cpp)
TARGET_LINK_LIBRARIES(complexp_test ${MANDATORY_LIBRARIES} piranha_static)

ADD_EXECUTABLE(tass_r6_test tass_r6_test.cpp)
TARGET_LINK_LIBRARIES(tass_r6_test ${MANDATORY_LIBRARIES} piranha_static)
