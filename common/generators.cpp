#include "generators.h"

namespace LIB_DTS2_NAMESPACE {

std::vector< std::pair< uint32_t, uint32_t > > generateUIntPairs(std::size_t count, const std::pair< uint32_t, uint32_t >& minMaxFirst, const std::pair< uint32_t, uint32_t >& minMaxSecond) {
	std::vector< std::pair< uint32_t, uint32_t > > tmp;
	generateUIntPairs(count, minMaxFirst, minMaxSecond, std::back_inserter(tmp));
	return tmp;
}

std::vector< std::pair< uint32_t, uint32_t > > generateUIntPairs(std::size_t count, const std::pair< uint32_t, uint32_t >& minMax) {
	return generateUIntPairs(count, minMax, minMax);
}


} //end namespace dts2