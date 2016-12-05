cmake_minimum_required(VERSION 2.8)
project(dts2common)

include_directories(
	${PROJECT_GLOBAL_SOURCE_DIR}
	${PROJECT_GLOBAL_SOURCE_DIR}/include
	${CMAKE_CURRENT_BINARY_DIR}
	${CMAKE_CURRENT_SOURCE_DIR}
	${CMAKE_SOURCE_DIR}
	${CMAKE_SOURCE_DIR}/include
)

set(COMMON_LIB_SOURCES_CPP
	generators.cpp
)

add_library(${PROJECT_NAME} STATIC ${COMMON_LIB_SOURCES_CPP})
add_dependencies(${PROJECT_NAME} dts2)