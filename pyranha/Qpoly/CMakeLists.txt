SET(PYRANA_QPOLY_HEADERS 
    ../../src/manipulators/qpoly.h
    ../series_instantiations.h
    ../exceptions.h
)

ADD_LIBRARY(_Qpoly SHARED Qpoly.cpp ${PYRANA_QPOLY_HEADERS})
PYRANHA_MANIPULATOR(_Qpoly)
