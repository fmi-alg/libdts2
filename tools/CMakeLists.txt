cmake_minimum_required(VERSION 3.9)
project(dts2tools)

set(TOOLS_LIB_SOURCES_CPP
	types.cpp
)

set(LINK_LIBS
	dts2
)

add_library(${PROJECT_NAME} STATIC ${TOOLS_LIB_SOURCES_CPP})

MACRO(ADD_TOOLS_TARGET _name _sourcefiles)
	add_executable("${PROJECT_NAME}_${_name}" "${_sourcefiles}")
	target_link_libraries("${PROJECT_NAME}_${_name}" ${LINK_LIBS})
	add_dependencies("${PROJECT_NAME}_${_name}" ratsstools)
	set_target_properties("${PROJECT_NAME}_${_name}" PROPERTIES OUTPUT_NAME ${_name})
ENDMACRO(ADD_TOOLS_TARGET)

ADD_TOOLS_TARGET(triang triang.cpp)
