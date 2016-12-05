#ifndef LIB_DTS2_COMMON_GENERATORS_H
#define LIB_DTS2_COMMON_GENERATORS_H
#pragma once
#include <libdts2/constants.h>

#include <vector>
#include <random>

namespace LIB_DTS2_NAMESPACE {

template<typename T_OUTPUT_ITERATOR>
void generateUIntPairs(std::size_t count, const std::pair<uint32_t, uint32_t> & minMaxFirst, const std::pair<uint32_t, uint32_t> & minMaxSecond, T_OUTPUT_ITERATOR out);

std::vector< std::pair<uint32_t, uint32_t> > generateUIntPairs(std::size_t count, const std::pair<uint32_t, uint32_t> & minMaxFirst, const std::pair<uint32_t, uint32_t> & minMaxSecond);
std::vector< std::pair<uint32_t, uint32_t> > generateUIntPairs(std::size_t count, const std::pair<uint32_t, uint32_t> & minMax);

}

//definitions
namespace LIB_DTS2_NAMESPACE {

template<typename T_OUTPUT_ITERATOR>
void generateUIntPairs(std::size_t count, const std::pair< uint32_t, uint32_t >& minMaxFirst, const std::pair< uint32_t, uint32_t >& minMaxSecond, T_OUTPUT_ITERATOR out) {
	std::uniform_int_distribution<uint32_t> dsF(minMaxFirst.first, minMaxFirst.second);
	std::uniform_int_distribution<uint32_t> dsS(minMaxSecond.first, minMaxSecond.second);
	std::default_random_engine gen;
	
	std::pair<uint32_t, uint32_t> tmp;
	for(std::size_t i(0); i < count; ++i) {
		tmp.first = dsF(gen);
		tmp.second = dsS(gen);
		*out = tmp;
		++out;
	}
}

}

#endif