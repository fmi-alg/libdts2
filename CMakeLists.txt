cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH
	${CMAKE_SOURCE_DIR}/cmake
	${CMAKE_CURRENT_SOURCE_DIR}/cmake
)

project(dts2)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++0x")

find_package(CGAL REQUIRED Core)
add_subdirectory(vendor/libratss libratss)

include_directories(
	${CMAKE_SOURCE_DIR}
	${CMAKE_CURRENT_SOURCE_DIR}
	${CMAKE_SOURCE_DIR}/include
	${CMAKE_CURRENT_SOURCE_DIR}/include
	${LIBRATSS_INCLUDE_DIR}
)

message("LIBRATS_INCLUDE_DIR=${LIBRATSS_INCLUDE_DIR}")

get_directory_property(hasParent PARENT_DIRECTORY)
#export correct link libs and include dirs
if(hasParent)
	SET(LIBDTS2_INCLUDE_DIR
		"${CMAKE_CURRENT_SOURCE_DIR}/include"
		${LIBRATSS_INCLUDE_DIR}
		PARENT_SCOPE)
	
	SET(LIBDTS2_LIBRARIES
		"${CMAKE_CURRENT_BINARY_DIR}/libdts2.a"
		PARENT_SCOPE)
	
	SET(LIBDTS2_LINK_LIBRARIES
		"${CMAKE_CURRENT_BINARY_DIR}/libdts2.a"
		${LIBRATSS_LIBRARIES}
		"CGAL"
		"CGAL_Core"
		"gmpxx"
		"mpfr"
		"gmp"
		"boost_thread"
		"boost_system"
		"pthread"
		PARENT_SCOPE
	)
endif()#hasParent

set(LIB_SOURCES_CPP
	src/Constrained_delaunay_triangulation_s2.cpp
)

add_library(${PROJECT_NAME} STATIC
	${LIB_SOURCES_CPP}
)
add_dependencies(${PROJECT_NAME} ratss)

add_subdirectory(common ${PROJECT_NAME}common EXCLUDE_FROM_ALL)
add_subdirectory(bench ${PROJECT_NAME}bench EXCLUDE_FROM_ALL)
add_subdirectory(tests ${PROJECT_NAME}tests EXCLUDE_FROM_ALL)

