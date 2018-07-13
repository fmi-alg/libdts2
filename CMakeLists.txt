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

cmake_minimum_required(VERSION 3.9)
project(dts2)

set(CMAKE_MODULE_PATH
	"${CMAKE_CURRENT_SOURCE_DIR}/cmake"
	"${CMAKE_CURRENT_SOURCE_DIR}/vendor/libratss/cmake"
)

set(MY_CXX_FLAGS "-Wall -Wextra -Wstrict-overflow=3")
set(MY_C_FLAGS "-Wall -Wextra -Wstrict-overflow=3")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MY_CXX_FLAGS}")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MY_C_FLAGS}")

set(MY_INCLUDE_DIRS
	"${CMAKE_CURRENT_SOURCE_DIR}/include"
)

find_package(CGAL REQUIRED Core)
add_subdirectory(vendor/libratss libratss)

set(LIB_SOURCES_CPP
	src/Constrained_delaunay_triangulation_s2.cpp
)

add_library(${PROJECT_NAME} STATIC
	${LIB_SOURCES_CPP}
)
target_link_libraries(${PROJECT_NAME} ratss)
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_14)
target_include_directories(${PROJECT_NAME} PUBLIC ${MY_INCLUDE_DIRS})


add_subdirectory(common ${PROJECT_NAME}common EXCLUDE_FROM_ALL)
add_subdirectory(bench ${PROJECT_NAME}bench EXCLUDE_FROM_ALL)
add_subdirectory(tests ${PROJECT_NAME}tests EXCLUDE_FROM_ALL)
add_subdirectory(tools ${PROJECT_NAME}tools EXCLUDE_FROM_ALL)
add_subdirectory(examples ${PROJECT_NAME}examples EXCLUDE_FROM_ALL)
