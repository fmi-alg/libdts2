cmake_minimum_required(VERSION 3.5)
project(dts2examples)

include_directories(
	${LIBDTS2_INCLUDE_DIR}
)

message(STATUS "dts2examples include dir: ${LIBDTS2_INCLUDE_DIR}")

set(LINK_LIBS
	${LIBDTS2_LINK_LIBRARIES}
)

MACRO(ADD_EXAMPLES_TARGET _name _sourcefiles)
	add_executable("${PROJECT_NAME}_${_name}" "${_sourcefiles}")
	target_link_libraries("${PROJECT_NAME}_${_name}" ${LINK_LIBS})
	add_dependencies("${PROJECT_NAME}_${_name}" dts2 ratss ratsscommon)
	set_target_properties("${PROJECT_NAME}_${_name}" PROPERTIES OUTPUT_NAME ${_name})
ENDMACRO(ADD_EXAMPLES_TARGET)

ADD_EXAMPLES_TARGET(simple simple.cpp)
