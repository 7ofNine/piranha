SET(PYRANHA_QQPOLY_HEADERS
    ../../src/manipulators/qqpoly.h
    ../series_instantiations.h
    ../exceptions.h
)

ADD_LIBRARY(_Qqpoly SHARED Qqpoly.cpp ${PYRANHA_QQPOLY_HEADERS})
PYRANHA_MANIPULATOR(_Qqpoly)
