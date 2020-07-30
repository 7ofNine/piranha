SET(PYRANHA_QQPS_HEADERS
    ../../src/manipulators/qqps.h
    ../series_instantiations.h
    ../exceptions.h
)

SET(PYRANHA_QQPS_PYTHON
    __init__.py
    )
SOURCE_GROUP(TREE ${CMAKE_CURRENT_LIST_DIR} PREFIX "Python" FILES ${PYRANHA_QQPS_PYTHON})
ADD_LIBRARY(_Qqps SHARED Qqps.cpp ${PYRANHA_QQPS_HEADERS} ${PYRANHA_QQPS_PYTHON})
PYRANHA_MANIPULATOR(_Qqps)
