cmake_minimum_required(VERSION 3.9)
project(dts2tools)

find_package(OpenMP)

set(TOOLS_LIB_SOURCES_CPP
	types.cpp
)

add_library(${PROJECT_NAME} STATIC ${TOOLS_LIB_SOURCES_CPP})
target_link_libraries(${PROJECT_NAME} dts2)

MACRO(ADD_TOOLS_TARGET _name _sourcefiles)
	add_executable("${PROJECT_NAME}_${_name}" "${_sourcefiles}")
	target_link_libraries("${PROJECT_NAME}_${_name}" dts2tools)
	set_target_properties("${PROJECT_NAME}_${_name}" PROPERTIES OUTPUT_NAME ${_name})
	if (OpenMP_CXX_FOUND)
		target_link_libraries("${PROJECT_NAME}_${_name}" OpenMP::OpenMP_CXX)
	endif()
ENDMACRO(ADD_TOOLS_TARGET)

ADD_TOOLS_TARGET(triang triang.cpp)
