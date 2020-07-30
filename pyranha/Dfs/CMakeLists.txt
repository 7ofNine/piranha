SET(PYRANHA_DFS_HEADERS 
	../../src/manipulators/dfs.h
	../series_instantiations.h
	../exceptions.h
	)
SET(PYRANHA_DFS_PYTHON
    __init__.py
    )
SOURCE_GROUP(TREE ${CMAKE_CURRENT_LIST_DIR} PREFIX "Python" FILES ${PYRANHA_DFS_PYTHON})    
ADD_LIBRARY(_Dfs SHARED Dfs.cpp ${PYRANHA_DFS_HEADERS} ${PYRANHA_DFS_PYTHON} )
PYRANHA_MANIPULATOR(_Dfs)
