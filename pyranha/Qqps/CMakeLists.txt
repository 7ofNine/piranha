SET(PYRANHA_QQPS_HEADERS
    ../../src/manipulators/qqps.h
    ../series_instantiations.h
    ../exceptions.h
)

ADD_LIBRARY(_Qqps SHARED Qqps.cpp ${PYRANHA_QQPS_HEADERS})
PYRANHA_MANIPULATOR(_Qqps)
