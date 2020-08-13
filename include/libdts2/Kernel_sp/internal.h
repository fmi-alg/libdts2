#pragma once
#ifndef LIB_DTS2_KERNEL_SP_INTERNAL_H
#define LIB_DTS2_KERNEL_SP_INTERNAL_H

#include <boost/multiprecision/cpp_int.hpp>

#include <libdts2/vendor/ttmath/ttmathint.h>
// #define LIB_DTS_2_KERNEL_SP_USE_TTMATH
namespace LIB_DTS2_NAMESPACE::detail::Kernel_sp {

#ifdef LIB_DTS_2_KERNEL_SP_USE_TTMATH

template<int T_BITS>
struct IntegerTypeFromBits {
	static constexpr int words = std::max<int>(1, (T_BITS/64) + int(T_BITS%64!=0));
	using type = ttmath::Int<words>;
};

template<>
struct IntegerTypeFromBits<32> {
	using type = int64_t;
};

template<>
struct IntegerTypeFromBits<64> {
	using type = int64_t;
};

#else
	
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

#endif

template<int T_BITS>
struct AlignedIntegerTypeFromBits {
	static constexpr int requested_bits = T_BITS;
	static constexpr int bits = (requested_bits/64 + int((requested_bits%64)!=0))*64;
	using type = typename IntegerTypeFromBits<bits>::type;
};

template<typename T_RT, int T_MAX_BITS = std::numeric_limits<T_RT>::digits>
struct Numerator_Denominator_1_3 {
	using RT = T_RT;
	static constexpr int max_bits = T_MAX_BITS;
	RT num;
	RT den;
	Numerator_Denominator_1_3() {}
	Numerator_Denominator_1_3(RT const & num, RT const & den) : num(num), den(den) {}
	Numerator_Denominator_1_3(Numerator_Denominator_1_3 const &) = default;
	Numerator_Denominator_1_3(Numerator_Denominator_1_3 &&) = default;
	Numerator_Denominator_1_3 & operator=(Numerator_Denominator_1_3 const &) = default;
	Numerator_Denominator_1_3 & operator=(Numerator_Denominator_1_3 &&) = default;
};

template<typename T_RT, int T_MAX_BITS = std::numeric_limits<T_RT>::digits>
struct Numerators_3 {
	using RT = T_RT;
	static constexpr int max_bits = T_MAX_BITS;
	RT x;
	RT y;
	RT z;
	Numerators_3() {}
	Numerators_3(RT const & x, RT const & y, RT const & z) : x(x), y(y), z(z) {}
	Numerators_3(Numerators_3 const &) = default;
	Numerators_3(Numerators_3 &&) = default;
	Numerators_3 & operator=(Numerators_3 const &) = default;
	Numerators_3 & operator=(Numerators_3 &&) = default;
};

} // LIB_DTS2_NAMESPACE::detail::Kernel_sp

#endif
