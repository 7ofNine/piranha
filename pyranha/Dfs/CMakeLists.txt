SET(DFS_HEADERS ../../src/manipulators/dfs.h ../series_instantiations.h ../exceptions.h)

ADD_LIBRARY(_Dfs SHARED Dfs.cpp ${GUT_PYRANHA_ADDED_SRC} ${DFS_HEADERS})
PYRANHA_MANIPULATOR(_Dfs)
