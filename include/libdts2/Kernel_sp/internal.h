#pragma once
#ifndef LIB_DTS2_KERNEL_SP_INTERNAL_H
#define LIB_DTS2_KERNEL_SP_INTERNAL_H

#include <boost/multiprecision/cpp_int.hpp>

namespace LIB_DTS2_NAMESPACE::detail::Kernel_sp {

template<int T_BITS>
struct IntegerTypeFromBits {
	using type = boost::multiprecision::number<
		boost::multiprecision::cpp_int_backend<
			T_BITS,
			T_BITS,
			boost::multiprecision::signed_magnitude,
			boost::multiprecision::unchecked,
			void
		>
	>;
};

template<>
struct IntegerTypeFromBits<32> {
	using type = int32_t;
};

template<>
struct IntegerTypeFromBits<64> {
	using type = int64_t;
};

template<>
struct IntegerTypeFromBits<128> {
	using type = __int128_t;
};

template<int T_BITS>
struct AlignedIntegerTypeFromBits {
	static constexpr int requested_bits = T_BITS;
	static constexpr int bits = (requested_bits/64 + int((requested_bits%64)!=0))*64;
	using type = typename IntegerTypeFromBits<bits>::type;
};

template<typename T_RT>
struct Numerators_3 {
	using RT = T_RT;
	RT x;
	RT y;
	RT z;
	Numerators_3() {}
	Numerators_3(RT const & x, RT const & y, RT const & z) : x(x), y(y), z(z) {}
	Numerators_3(Numerators_3 const &) = default;
	Numerators_3(Numerators_3 &&) = default;
};

} // LIB_DTS2_NAMESPACE::detail::Kernel_sp

#endif
