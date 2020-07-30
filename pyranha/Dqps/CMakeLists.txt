SET(PYRANHA_DQPS_HEADERS 
    ../../src/manipulators/dqps.h
    ../series_instantiations.h
    ../exceptions.h
)
SET(PYRANHA_DQPS_PYTHON
    __init__.py
    )
SOURCE_GROUP(TREE ${CMAKE_CURRENT_LIST_DIR} PREFIX "Python" FILES ${PYRANHA_DQPS_PYTHON}) 
ADD_LIBRARY(_Dqps SHARED Dqps.cpp ${PYRANHA_DQPS_HEADERS} ${PYRANHA_DQPS_PYTHON})
PYRANHA_MANIPULATOR(_Dqps)
