cmake_minimum_required(VERSION 2.8)
project(dts2tools)

include_directories(
	${LIBDTS2_INCLUDE_DIR}
)

set(TOOLS_LIB_SOURCES_CPP
	types.cpp
)

add_library(${PROJECT_NAME} STATIC ${TOOLS_LIB_SOURCES_CPP})
add_dependencies(${PROJECT_NAME} ratss ratsscommon)

set(LINK_LIBS
	${LIBRATSS_LINK_LIBRARIES}
)

MACRO(ADD_TOOLS_TARGET _name _sourcefiles)
	add_executable("${PROJECT_NAME}_${_name}" "${_sourcefiles}")
	target_link_libraries("${PROJECT_NAME}_${_name}" ${LINK_LIBS})
	add_dependencies("${PROJECT_NAME}_${_name}" ratsstools)
	set_target_properties("${PROJECT_NAME}_${_name}" PROPERTIES OUTPUT_NAME ${_name})
ENDMACRO(ADD_TOOLS_TARGET)

ADD_TOOLS_TARGET(triang triang.cpp)