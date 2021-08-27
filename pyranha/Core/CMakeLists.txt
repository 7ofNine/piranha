SET(PYRANHA_CORE_HEADERS
    ../../src/core/base_classes/named_series_def.h
    ../../src/core/config.h
    ../../src/core/mp.h
    ../../src/core/Psym.h
    ../../src/core/settings.h
    ../../src/core/stats.h
    ../args_tuple.h
    ../boost_python_container_conversions.h
    ../commons.h
    ../exceptions.h
    ../mp_classes.h
)

SET(PYRANHA_CORE_PYTHON
    __init__.py
    impl.py
)
SOURCE_GROUP(TREE ${CMAKE_CURRENT_LIST_DIR} PREFIX "Python" FILES ${PYRANHA_CORE_PYTHON})
ADD_LIBRARY(_Core SHARED Core.cpp ${PYRANHA_CORE_HEADERS} ${PYRANHA_CORE_PYTHON})
PYTHON_CPPMODULE(_Core)
INSTALL(FILES __init__.py impl.py DESTINATION ${PYRANHA_INSTALL_PATH}/Core)
