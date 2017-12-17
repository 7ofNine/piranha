SET(PYRANHA_QPS_HEADERS 
    ../../src/manipulators/qps.h
    ../series_instantiations.h
    ../exceptions.h
)

ADD_LIBRARY(_Qps SHARED Qps.cpp ${PYRANHA_QPS_HEADERS})
PYRANHA_MANIPULATOR(_Qps)
