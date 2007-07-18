SET(SRC_LIST
  converter/arg_to_python_base.cpp
  converter/builtin_converters.cpp
  converter/from_python.cpp
  converter/registry.cpp
  converter/type_id.cpp
  object/class.cpp
  object/enum.cpp
  object/function.cpp
  object/inheritance.cpp
  object/iterator.cpp
  object/life_support.cpp
  object/pickle_support.cpp
  object/stl_iterator.cpp
  dict.cpp
  errors.cpp
  exec.cpp
  import.cpp
  list.cpp
  long.cpp
  module.cpp
  numeric.cpp
  object_operators.cpp
  object_protocol.cpp
  slice.cpp
  str.cpp
  tuple.cpp
  wrapper.cpp
  )

ADD_LIBRARY(${BOOST_PYTHON_LIBRARY} SHARED ${SRC_LIST})

TARGET_LINK_LIBRARIES(${BOOST_PYTHON_LIBRARY} ${MANDATORY_LIBRARIES})
INSTALL(TARGETS ${BOOST_PYTHON_LIBRARY}
  RUNTIME DESTINATION ${LIB_INSTALL_PATH}
  LIBRARY DESTINATION ${LIB_INSTALL_PATH}
)