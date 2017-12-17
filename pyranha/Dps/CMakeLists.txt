SET(PYRANHA_DPS_HEADERS 
    ../../src/manipulators/dfs.h
    ../../src/manipulators/dps.h
    ../series_instantiations.h
    ../exceptions.h
)
ADD_LIBRARY(_Dps SHARED Dps.cpp ${PYRANHA_DPS_HEADERS})
PYRANHA_MANIPULATOR(_Dps)
