SET(PYRANHA_MATH_HEADERS 
    ../../src/core/math.h
    ../../src/core/mp.h
    ../exceptions.h
)

ADD_LIBRARY(_Math SHARED Math.cpp ${PYRANHA_MATH_HEADERS})
PYTHON_CPPMODULE(_Math)
INSTALL(FILES __init__.py DESTINATION ${PYRANHA_INSTALL_PATH}/Math/)
