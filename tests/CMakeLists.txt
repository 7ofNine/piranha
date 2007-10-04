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

ADD_EXECUTABLE(tass_r6_test tass_r6_test.cpp)
TARGET_LINK_LIBRARIES(tass_r6_test ${MANDATORY_LIBRARIES} piranha_static)

IF(ENABLE_SSE2)
  ADD_EXECUTABLE(gastest_elp3_simd gastest_elp3_simd.cpp)
  TARGET_LINK_LIBRARIES(gastest_elp3_simd ${MANDATORY_LIBRARIES} piranha_static)
ENDIF(ENABLE_SSE2)
