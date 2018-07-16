cmake_minimum_required(VERSION 3.9)
project(dts2bench)

MACRO(ADD_BENCH_TARGET _name _sourcefiles)
	add_executable("${_name}" "${_sourcefiles}")
	target_link_libraries("${_name}" dts2 dts2common ratsscommon)
ENDMACRO(ADD_BENCH_TARGET)
