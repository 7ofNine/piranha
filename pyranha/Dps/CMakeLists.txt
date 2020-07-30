SET(PYRANHA_DPS_HEADERS 
    ../../src/manipulators/dfs.h
    ../../src/manipulators/dps.h
    ../series_instantiations.h
    ../exceptions.h
)
SET(PYRANHA_DPS_PYTHON
    __init__.py
    )
SOURCE_GROUP(TREE ${CMAKE_CURRENT_LIST_DIR} PREFIX "Python" FILES ${PYRANHA_DPS_PYTHON}) 
ADD_LIBRARY(_Dps SHARED Dps.cpp ${PYRANHA_DPS_HEADERS} ${PYRANHA_DPS_PYTHON})
PYRANHA_MANIPULATOR(_Dps)
