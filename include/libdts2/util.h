#pragma once
#ifndef LIB_DTS2_UTIL_H
#define LIB_DTS2_UTIL_H

#include <libdts2/constants.h>

namespace LIB_DTS2_NAMESPACE {
namespace detail {
	
///Identity pass through cast
///This makes sure that casting the same class to the same class does not incur any overhead
template<typename T_TARGET_TYPE, typename T_SOURCE_TYPE>
struct ipt_static_cast_impl {
	using SourceType = T_SOURCE_TYPE;
	using TargetType = T_TARGET_TYPE;
	static TargetType c(SourceType const & s) { return static_cast<TargetType>(s); }
	static TargetType c(SourceType & s) { return static_cast<TargetType>(s); }
};

template<typename T_TARGET_TYPE>
struct ipt_static_cast_impl<T_TARGET_TYPE, T_TARGET_TYPE> {
	using SourceType = T_TARGET_TYPE;
	using TargetType = T_TARGET_TYPE;
	static TargetType const & c(SourceType const & s) { return s; }
	static TargetType & c(SourceType & s) { return s; }
};

template<typename T_TARGET_TYPE, typename T_SOURCE_TYPE>
auto ipt_static_cast(T_SOURCE_TYPE const & v) {
	return ipt_static_cast_impl<T_TARGET_TYPE, T_SOURCE_TYPE>::c(v);
}

template<typename T_TARGET_TYPE, typename T_SOURCE_TYPE>
auto ipt_static_cast(T_SOURCE_TYPE & v) {
	return ipt_static_cast_impl<T_TARGET_TYPE, T_SOURCE_TYPE>::c(v);
}
	
}//end namespace detail
	
}//end namespace dts2

#endif
