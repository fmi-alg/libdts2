cmake_minimum_required(VERSION 3.9)
project(dts2common)

set(COMMON_LIB_SOURCES_CPP
	generators.cpp
)

add_library(${PROJECT_NAME} STATIC ${COMMON_LIB_SOURCES_CPP})
target_link_libraries(${PROJECT_NAME} dts2)
