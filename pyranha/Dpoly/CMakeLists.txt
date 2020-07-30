SET(PYRANHA_DPOLY_HEADERS
    ../../src/manipulators/dpoly.h
    ../series_instantiations.h
    ../exceptions.h
)
SET(PYRANHA_DPOLY_PYTHON
    __init__.py
    )
SOURCE_GROUP(TREE ${CMAKE_CURRENT_LIST_DIR} PREFIX "Python" FILES ${PYRANHA_DPOLY_PYTHON}) 
ADD_LIBRARY(_Dpoly SHARED Dpoly.cpp ${PYRANHA_DPOLY_HEADERS} ${PYRANHA_DPOLY_PYTHON})
PYRANHA_MANIPULATOR(_Dpoly)
