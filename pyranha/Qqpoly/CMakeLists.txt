SET(PYRANHA_QQPOLY_HEADERS
    ../../src/manipulators/qqpoly.h
    ../series_instantiations.h
    ../exceptions.h
)

SET(PYRANHA_QQPOLY_PYTHON
    __init__.py
    )
SOURCE_GROUP(TREE ${CMAKE_CURRENT_LIST_DIR} PREFIX "Python" FILES ${PYRANHA_QQPOLY_PYTHON})
ADD_LIBRARY(_Qqpoly SHARED Qqpoly.cpp ${PYRANHA_QQPOLY_HEADERS} ${PYRANHA_QQPOLY_PYTHON})
PYRANHA_MANIPULATOR(_Qqpoly)
