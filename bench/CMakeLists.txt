cmake_minimum_required(VERSION 2.8)
project(dts2bench)

include_directories(
	${PROJECT_GLOBAL_SOURCE_DIR}
	${CMAKE_CURRENT_BINARY_DIR}
	${CMAKE_CURRENT_SOURCE_DIR}
	${CMAKE_SOURCE_DIR}
)

set(LINK_LIBS
	dts2common
	ratsscommon
	${LIBDTS2_LINK_LIBRARIES}
	${CPPUNIT_LIBRARY}
)

# add_dependencies(${PROJECT_NAME} dts2 common)

MACRO(ADD_BENCH_TARGET _name _sourcefiles)
	add_executable("${_name}" "${_sourcefiles}")
	target_link_libraries("${_name}" ${LINK_LIBS})
	add_dependencies("${_name}" dts2 common)
ENDMACRO(ADD_BENCH_TARGET)
