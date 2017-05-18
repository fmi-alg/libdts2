cmake_minimum_required(VERSION 2.8)
project(dts2tests)
find_package(CPPUNIT REQUIRED)

include_directories(
	${PROJECT_GLOBAL_SOURCE_DIR}
	${PROJECT_GLOBAL_SOURCE_DIR}/include
	${CMAKE_CURRENT_BINARY_DIR}
	${CMAKE_CURRENT_SOURCE_DIR}
	${CMAKE_SOURCE_DIR}
	${CMAKE_SOURCE_DIR}/include
)

set(TEST_LIB_SOURCES_CPP
	TestBase.cpp
)

add_library(${PROJECT_NAME} STATIC ${TEST_LIB_SOURCES_CPP})
add_dependencies(${PROJECT_NAME} dts2 ratsscommon)

MACRO(ADD_TEST_TARGET_SINGLE _name)
	add_executable("${PROJECT_NAME}_${_name}" "${_name}.cpp")
	target_link_libraries("${PROJECT_NAME}_${_name}" ${TEST_LIB_LINK_LIBS})
	add_dependencies("${PROJECT_NAME}_${_name}" dts2tests)
	set_target_properties("${PROJECT_NAME}_${_name}" PROPERTIES OUTPUT_NAME ${_name})
ENDMACRO(ADD_TEST_TARGET_SINGLE)

set(TEST_LIB_LINK_LIBS
	dts2tests
	ratsscommon
	${LIBDTS2_LINK_LIBRARIES}
	${CPPUNIT_LIBRARY}
)

ADD_TEST_TARGET_SINGLE(dt)
ADD_TEST_TARGET_SINGLE(cdt)