cmake_minimum_required(VERSION 3.9)
project(dts2examples)

MACRO(ADD_EXAMPLES_TARGET _name _sourcefiles)
	add_executable("${PROJECT_NAME}_${_name}" "${_sourcefiles}")
	target_link_libraries("${PROJECT_NAME}_${_name}" dts2)
	set_target_properties("${PROJECT_NAME}_${_name}" PROPERTIES OUTPUT_NAME ${_name})
ENDMACRO(ADD_EXAMPLES_TARGET)

ADD_EXAMPLES_TARGET(simple simple.cpp)
