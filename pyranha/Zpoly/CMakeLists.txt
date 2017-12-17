SET(PYRANHA_ZPOLY_HEADERS 
    ../../src/manipulators/zpoly.h
    ../series_instantiations.h
    ../exceptions.h
)

ADD_LIBRARY(_Zpoly SHARED Zpoly.cpp ${PYRANHA_ZPOLY_HEADERS})
PYRANHA_MANIPULATOR(_Zpoly)
