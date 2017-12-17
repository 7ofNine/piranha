SET(PYRANHA_DQPS_HEADERS 
    ../../src/manipulators/dqps.h
    ../series_instantiations.h
    ../exceptions.h
)

ADD_LIBRARY(_Dqps SHARED Dqps.cpp ${PYRANHA_DQPS_HEADERS})
PYRANHA_MANIPULATOR(_Dqps)
