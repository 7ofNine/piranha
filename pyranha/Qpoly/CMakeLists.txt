SET(PYRANA_QPOLY_HEADERS 
    ../../src/manipulators/qpoly.h
    ../series_instantiations.h
    ../exceptions.h
)

SET(PYRANHA_QPOLY_PYTHON
    __init__.py
    )
SOURCE_GROUP(TREE ${CMAKE_CURRENT_LIST_DIR} PREFIX "Python" FILES ${PYRANHA_QPOLY_PYTHON})
ADD_LIBRARY(_Qpoly SHARED Qpoly.cpp ${PYRANA_QPOLY_HEADERS} ${PYRANHA_QPOLY_PYTHON})
PYRANHA_MANIPULATOR(_Qpoly)
