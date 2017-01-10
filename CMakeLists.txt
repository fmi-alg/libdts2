# This is the build file to build libdts2
# It can be build as a stand-alone project as well as a sub-project
# If this is intended as a sub-project then make sure to add the appropriate compile definitions to your project
# This build file exports the following variables into its parent scope:
# LIBDTS2_INCLUDE_DIR - The include directory
# LIBDTS2_LIBRARIES - libraries that need to be linker in order for libdts2 to work
# LIBDTS2_COMPILE_DEFINITIONS - definitions that are needed
# 
# You may simply do the following:
# 
# include_directories(${LIBDTS2_INCLUDE_DIR})
# 
# add_executable(${PROJECT_NAME} ${SOURCES_CPP}
# target_link_libraries(${PROJECT_NAME} ${LIBDTS2_LIBRARIES})
# target_compile_definitions(${PROJECT_NAME} PUBLIC ${LIBDTS2_COMPILE_DEFINITIONS})
# add_dependencies(${PROJECT_NAME} dts2)

cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH
	${CMAKE_SOURCE_DIR}/cmake
	${CMAKE_CURRENT_SOURCE_DIR}/cmake
)

project(dts2)

set(MY_CXX_FLAGS "-std=gnu++0x")
set(MY_C_FLAGS "-std=c99")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MY_CXX_FLAGS}")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MY_C_FLAGS}")

set(MY_LTO_CXX_FLAGS "${MY_CXX_FLAGS} -DNDEBUG -O3 -march=native -flto -ffat-lto-objects -frounding-math")
set(MY_LTO_C_FLAGS "${MY_C_FLAGS} -DNDEBUG -O3 -march=native -flto -ffat-lto-objects -frounding-math")
set(MY_LTO_LINKER_FLAGS "")

SET(CMAKE_CXX_FLAGS_LTO
	"${MY_LTO_CXX_FLAGS}"
	CACHE STRING "Flags used by the C++ compiler during lto builds."
	FORCE
)
SET(CMAKE_C_FLAGS_LTO
	"${MY_PROFILE_GENERATE_C_FLAGS}"
	CACHE STRING "Flags used by the C compiler during lto builds."
	FORCE
)
SET(CMAKE_EXE_LINKER_FLAGS_LTO
	"${MY_PROFILE_GENERATE_LINKER_FLAGS}"
	CACHE STRING "Flags used for linking binaries during lto builds."
	FORCE
)
SET(CMAKE_SHARED_LINKER_FLAGS_LTO
	""
	CACHE STRING "Flags used by the shared libraries linker during lto builds."
	FORCE
)
MARK_AS_ADVANCED(
	CMAKE_CXX_FLAGS_LTO
	CMAKE_C_FLAGS_LTO
	CMAKE_EXE_LINKER_FLAGS_LTO
	CMAKE_SHARED_LINKER_FLAGS_LTO
)


find_package(CGAL REQUIRED Core)
add_subdirectory(vendor/libratss libratss)

include_directories(
# 	${CMAKE_SOURCE_DIR}
# 	${CMAKE_CURRENT_SOURCE_DIR}
# 	${CMAKE_SOURCE_DIR}/include
	${CMAKE_CURRENT_SOURCE_DIR}/include
	${LIBRATSS_INCLUDE_DIR}
)

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

set(LIBDTS2_COMPILE_DEFINITIONS
	${LIBRATSS_COMPILE_DEFINITIONS}
)

get_directory_property(hasParent PARENT_DIRECTORY)
if(hasParent)
	set(LIBDTS2_INCLUDE_DIR
		${LIBDTS2_INCLUDE_DIR}
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
	
	set(LIBDTS2_COMPILE_DEFINITIONS
		${LIBDTS2_COMPILE_DEFINITIONS}
		PARENT_SCOPE
	)
endif(hasParent)

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

