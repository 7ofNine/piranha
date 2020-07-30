SET(PYRANHA_TRUNCATORS_HEADERS
    ../../src/core/truncators/degree.h
    ../../src/core/truncators/norm.h
    ../../src/core/truncators/truncators.h
    ../../src/core/mp.h
    ../commons.h
    ../exceptions.h
)

SET(PYRANHA_TRUNCATORS_PYTHON
    __init__.py
    detail.py
    )
SOURCE_GROUP(TREE ${CMAKE_CURRENT_LIST_DIR} PREFIX "Python" FILES ${PYRANHA_TRUNCATORS_PYTHON})
ADD_LIBRARY(_Truncators SHARED Truncators.cpp ${PYRANHA_TRUNCATORS_HEADERS} ${PYRANHA_TRUNCATORS_PYTHON})
PYTHON_CPPMODULE(_Truncators)
INSTALL(FILES __init__.py detail.py DESTINATION ${PYRANHA_INSTALL_PATH}/Truncators/)
