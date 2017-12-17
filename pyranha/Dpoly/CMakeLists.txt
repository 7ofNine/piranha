SET(PYRANHA_DPOLY_HEADERS
    ../../src/manipulators/dpoly.h
    ../series_instantiations.h
    ../exceptions.h
)
ADD_LIBRARY(_Dpoly SHARED Dpoly.cpp ${PYRANHA_DPOLY_HEADERS})
PYRANHA_MANIPULATOR(_Dpoly)
