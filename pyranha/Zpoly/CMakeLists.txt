SET(PYRANHA_ZPOLY_HEADERS 
    ../../src/manipulators/zpoly.h
    ../series_instantiations.h
    ../exceptions.h
)

SET(PYRANHA_ZPOLY_PYTHON
    __init__.py
    )
SOURCE_GROUP(TREE ${CMAKE_CURRENT_LIST_DIR} PREFIX "Python" FILES ${PYRANHA_ZPOLY_PYTHON})
ADD_LIBRARY(_Zpoly SHARED Zpoly.cpp ${PYRANHA_ZPOLY_HEADERS} ${PYRANHA_ZPOLY_PYTHON})
PYRANHA_MANIPULATOR(_Zpoly)
