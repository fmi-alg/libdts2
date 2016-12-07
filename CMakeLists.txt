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

message("LIBRATSS_COMPILE_DEFINITIONS=${LIBRATSS_COMPILE_DEFINITIONS}")

set(LIBDTS2_INCLUDE_DIR
	"${CMAKE_CURRENT_SOURCE_DIR}/include"
	${LIBRATSS_INCLUDE_DIR}
)

set(LIBDTS2_LIBRARIES
	dts2
)

set(LIBDTS2_LINK_LIBRARIES
	dts2
	${LIBRATSS_LINK_LIBRARIES}
)

get_directory_property(hasParent PARENT_DIRECTORY)
if(hasParent)
	set(LIBDTS2_INCLUDE_DIR
		${LIBDTS2_LINK_LIBRARIES}
		PARENT_SCOPE
	)
	
	set(LIBDTS2_LIBRARIES
		${LIBDTS2_LIBRARIES}
		PARENT_SCOPE
	)
	
	set(LIBDTS2_LINK_LIBRARIES
		${LIBDTS2_LINK_LIBRARIES}
		PARENT_SCOPE
	)
endif(hasParent)

set(LIB_SOURCES_CPP
	src/Constrained_delaunay_triangulation_s2.cpp
)

add_library(${PROJECT_NAME} STATIC
	${LIB_SOURCES_CPP}
)
target_compile_definitions(${PROJECT_NAME} PUBLIC ${LIBRATSS_COMPILE_DEFINITIONS})
add_dependencies(${PROJECT_NAME} ratss)

add_subdirectory(common ${PROJECT_NAME}common EXCLUDE_FROM_ALL)
add_subdirectory(bench ${PROJECT_NAME}bench EXCLUDE_FROM_ALL)
add_subdirectory(tests ${PROJECT_NAME}tests EXCLUDE_FROM_ALL)

