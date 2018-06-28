# This is the build file to build libdts2
# It can be build as a stand-alone project as well as a sub-project
# This build file exports the following variables into its parent scope:
# LIBDTS2_INCLUDE_DIR - The include directory
# LIBDTS2_LIBRARIES - libraries that need to be linker in order for libdts2 to work
# LIBDTS2_COMPILE_DEFINITIONS
# 
# You may simply do the following:
# 
# include_directories(${LIBDTS2_INCLUDE_DIR})
# 
# add_executable(${PROJECT_NAME} ${SOURCES_CPP}
# target_link_libraries(${PROJECT_NAME} ${LIBDTS2_LIBRARIES})
# add_dependencies(${PROJECT_NAME} dts2)

cmake_minimum_required(VERSION 3.5)
project(dts2)

set(CMAKE_MODULE_PATH
	"${CMAKE_CURRENT_SOURCE_DIR}/cmake"
	"${CMAKE_CURRENT_SOURCE_DIR}/vendor/libratss/cmake"
)

set(MY_CXX_FLAGS "-std=gnu++14 -Wall -Wextra -Wstrict-overflow=3")
set(MY_C_FLAGS "-std=c99 -Wall -Wextra -Wstrict-overflow=3")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MY_CXX_FLAGS}")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MY_C_FLAGS}")

find_package(CGAL REQUIRED Core)
add_subdirectory(vendor/libratss libratss)

set(LIBDTS2_INCLUDE_DIR
	"${CMAKE_CURRENT_SOURCE_DIR}/include"
	${LIBRATSS_INCLUDE_DIR}
	CACHE STRING "libdts2 include directories"
	FORCE
)

set(LIBDTS2_LIBRARIES
	${PROJECT_NAME}
	CACHE STRING "libdts2 libraries"
	FORCE
)

set(LIBDTS2_LINK_LIBRARIES
	${PROJECT_NAME}
	${LIBRATSS_LINK_LIBRARIES}
	CACHE STRING "libdts2 link libraries"
	FORCE
)

set(LIBDTS2_COMPILE_DEFINITIONS
	${LIBRATSS_COMPILE_DEFINITIONS}
	CACHE STRING "libdts2 compile definitions"
	FORCE
)

mark_as_advanced(
	LIBDTS2_INCLUDE_DIR
	LIBDTS2_LIBRARIES
	LIBDTS2_LINK_LIBRARIES
	LIBDTS2_COMPILE_DEFINITIONS
)

include_directories(
	${LIBDTS2_INCLUDE_DIR}
)

set(LIB_SOURCES_CPP
	src/Constrained_delaunay_triangulation_s2.cpp
)

add_library(${PROJECT_NAME} STATIC
	${LIB_SOURCES_CPP}
)

target_compile_definitions(${PROJECT_NAME} PUBLIC ${LIBDTS2_COMPILE_DEFINITIONS})
add_dependencies(${PROJECT_NAME} ratss)

add_subdirectory(common ${PROJECT_NAME}common EXCLUDE_FROM_ALL)
add_subdirectory(bench ${PROJECT_NAME}bench EXCLUDE_FROM_ALL)
add_subdirectory(tests ${PROJECT_NAME}tests EXCLUDE_FROM_ALL)
add_subdirectory(tools ${PROJECT_NAME}tools EXCLUDE_FROM_ALL)
