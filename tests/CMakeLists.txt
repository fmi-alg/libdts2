cmake_minimum_required(VERSION 3.9)
project(dts2tests)
find_package(CPPUNIT)

if (CPPUNIT_FOUND)

	set(TEST_LIB_SOURCES_CPP
		TestBase.cpp
	)

	add_library(${PROJECT_NAME} STATIC ${TEST_LIB_SOURCES_CPP})
	target_link_libraries(${PROJECT_NAME}
		ratsscommon
		dts2
		${CPPUNIT_LIBRARY}
	)

	MACRO(ADD_TEST_TARGET_SINGLE _name)
		add_executable("${PROJECT_NAME}_${_name}" "${_name}.cpp")
		target_link_libraries("${PROJECT_NAME}_${_name}" dts2tests)
		set_target_properties("${PROJECT_NAME}_${_name}" PROPERTIES OUTPUT_NAME ${_name})
	ENDMACRO(ADD_TEST_TARGET_SINGLE)

	ADD_TEST_TARGET_SINGLE(dt)
	ADD_TEST_TARGET_SINGLE(cdt)
	ADD_TEST_TARGET_SINGLE(kernel_sp)
else()
	message(WARNING "Unable to build tests due to missing cppunit")
endif()
