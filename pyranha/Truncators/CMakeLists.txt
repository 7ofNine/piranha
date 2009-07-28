ADD_LIBRARY(_Truncators SHARED Truncators.cpp)
PYTHON_CPPMODULE(_Truncators)
INSTALL(FILES __init__.py detail.py DESTINATION ${PYRANHA_INSTALL_PATH}/Truncators/)
