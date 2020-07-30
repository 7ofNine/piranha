SET(PYRANHA_QPS_HEADERS 
    ../../src/manipulators/qps.h
    ../series_instantiations.h
    ../exceptions.h
)

SET(PYRANHA_QPS_PYTHON
    __init__.py
    )
SOURCE_GROUP(TREE ${CMAKE_CURRENT_LIST_DIR} PREFIX "Python" FILES ${PYRANHA_QPS_PYTHON})
ADD_LIBRARY(_Qps SHARED Qps.cpp ${PYRANHA_QPS_HEADERS} ${PYRANHA_QPS_PYTHON})
PYRANHA_MANIPULATOR(_Qps)
