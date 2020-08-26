#pragma once
#ifndef LIB_DTS2_KERNEL_SP_INTERNAL_H
#define LIB_DTS2_KERNEL_SP_INTERNAL_H

#include <libdts2/constants.h>
#include <libdts2/util.h>
#include <libratss/Conversion.h>

// #define LIB_DTS_2_KERNEL_SP_USE_TTMATH
// #define LIB_DTS_2_KERNEL_SP_USE_MPZ_CLASS
// #define LIB_DTS_2_KERNEL_SP_USE_SIGNED_WIDER

#ifdef LIB_DTS_2_KERNEL_SP_USE_TTMATH
	#include <libdts2/vendor/ttmath/ttmathint.h>
	#include <libdts2/vendor/ttmath/ttmathuint.h>
#elif defined(LIB_DTS_2_KERNEL_SP_USE_MPZ_CLASS)
	#include <gmpxx.h>
#elif defined(LIB_DTS_2_KERNEL_SP_USE_SIGNED_WIDER)
	#include <libdts2/vendor/wider/signed_wider.h>
#else 
	#include <boost/multiprecision/cpp_int.hpp>
#endif

namespace LIB_DTS2_NAMESPACE::detail::Kernel_sp {
	
constexpr int align(int v, int alignment) { return (v/alignment) + int(v%alignment!=0); } 

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

//unsigned integer variants

template<int T_BITS>
struct UnsignedIntegerTypeFromBits {
	static constexpr int words = std::max<int>(1, align(T_BITS, 64));
	using type = ttmath::UInt<words>;
};

template<>
struct UnsignedIntegerTypeFromBits<32> {
	using type = uint64_t;
};

template<>
struct UnsignedIntegerTypeFromBits<64> {
	using type = uint64_t;
};

#elif defined(LIB_DTS_2_KERNEL_SP_USE_MPZ_CLASS)

template<int T_BITS>
struct IntegerTypeFromBits {
	using type = mpz_class;
};

template<>
struct IntegerTypeFromBits<32> {
	using type = int64_t;
};

template<>
struct IntegerTypeFromBits<64> {
	using type = int64_t;
};

//unsigned integer variants

template<int T_BITS>
struct UnsignedIntegerTypeFromBits {
	using type = mpz_class;
};

template<>
struct UnsignedIntegerTypeFromBits<32> {
	using type = uint64_t;
};

template<>
struct UnsignedIntegerTypeFromBits<64> {
	using type = uint64_t;
};

#elif defined(LIB_DTS_2_KERNEL_SP_USE_SIGNED_WIDER)

template<int T_BITS>
struct IntegerTypeFromBits;

template<>
struct IntegerTypeFromBits<32> {
	using type = int64_t;
};

template<>
struct IntegerTypeFromBits<64> {
	using type = int64_t;
};

template<>
struct IntegerTypeFromBits<128> {
	using type = SignedWider< Wider<uint64_t> >;
};

template<>
struct IntegerTypeFromBits<192> {
	using type = SignedWider< Wider< Wider<uint64_t> > >;
};

template<>
struct IntegerTypeFromBits<256> {
	using type = SignedWider< Wider< Wider<uint64_t> > >;
};

//unsigned integer variants

template<int T_BITS>
struct UnsignedIntegerTypeFromBits {
	using type = typename IntegerTypeFromBits<T_BITS>::type;
};

template<>
struct UnsignedIntegerTypeFromBits<32> {
	using type = uint64_t;
};

template<>
struct UnsignedIntegerTypeFromBits<64> {
	using type = uint64_t;
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

//unsigned integer variants

template<int T_BITS>
struct UnsignedIntegerTypeFromBits {
	using type = boost::multiprecision::number<
		boost::multiprecision::cpp_int_backend<
			T_BITS,
			T_BITS,
			boost::multiprecision::unsigned_magnitude,
			boost::multiprecision::unchecked,
			void
		>
	>;
};

template<>
struct UnsignedIntegerTypeFromBits<32> {
	using type = uint32_t;
};

template<>
struct UnsignedIntegerTypeFromBits<64> {
	using type = uint64_t;
};

template<>
struct UnsignedIntegerTypeFromBits<128> {
	using type = __uint128_t;
};

#endif

template<int T_DIGITS>
struct AlignedIntegerTypeFromDigits {
	static constexpr int requested_digits = T_DIGITS;
	static constexpr int requested_bits = requested_digits+1;
	static constexpr int bits = align(requested_bits, 64)*64;
	static constexpr int digits = bits-1;
	using type = typename IntegerTypeFromBits<bits>::type;
};

template<int T_DIGITS>
struct AlignedUnsignedIntegerTypeFromDigits {
	static constexpr int requested_digits = T_DIGITS;
	static constexpr int requested_bits = requested_digits;
	static constexpr int bits = align(requested_bits, 64)*64;
	static constexpr int digits = bits;
	using type = typename UnsignedIntegerTypeFromBits<bits>::type;
};

template<int T_DIGITS>
using AINT = typename AlignedIntegerTypeFromDigits<T_DIGITS>::type;

template<int T_DIGITS>
using AUINT = typename AlignedUnsignedIntegerTypeFromDigits<T_DIGITS>::type;

template<typename T_RT, int T_MAX_DIGITS = std::numeric_limits<T_RT>::digits>
struct Numerator_Denominator_1_3 {
	using RT = T_RT;
	static constexpr int max_digits = T_MAX_DIGITS;
	RT num;
	RT den;
	Numerator_Denominator_1_3() {}
	Numerator_Denominator_1_3(RT const & num, RT const & den) : num(num), den(den) {}
	Numerator_Denominator_1_3(Numerator_Denominator_1_3 const &) = default;
	Numerator_Denominator_1_3(Numerator_Denominator_1_3 &&) = default;
	Numerator_Denominator_1_3 & operator=(Numerator_Denominator_1_3 const &) = default;
	Numerator_Denominator_1_3 & operator=(Numerator_Denominator_1_3 &&) = default;
};

template<typename T_RT, int T_MAX_DIGITS = std::numeric_limits<T_RT>::digits>
struct Numerators_3 {
	using RT = T_RT;
	static constexpr int max_digits = T_MAX_DIGITS;
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

constexpr uint64_t log2(uint64_t v) {
	if (v == 0 || v == 1) {
		return 1;
	}
	else {
		return 1+log2(v/2);
	}
}

constexpr std::pair<const int, const bool> autoint_add_upper_bound(const int up1, const bool upl1, const int up2, const bool upl2) {
	using return_type = std::pair<const int, const bool>;
	if (upl1 && upl2) {
		return return_type(std::max(up1, up2)+1, true);
	}
	else if (!upl1 && !upl2) {
		if (int64_t(up1)+int64_t(up2) <= std::numeric_limits<int32_t>::max()) { //fits int
			return return_type(up1+up2, false);
		}
		else {
			return return_type(log2(uint64_t(up1)+uint64_t(up2)), true);
		}
	}
	else if (upl1 && !upl2) {
		if (up1 < 63) {
			const uint64_t up1_no_log = (static_cast<uint64_t>(1) << up1)-1;
			const uint64_t result_up_no_log = up1_no_log + uint64_t(up2);
			if (result_up_no_log < up1_no_log) { //overflow
				return return_type(up1+1, true);
			}
			else {
				if (result_up_no_log <= std::numeric_limits<int32_t>::max()) {
					return return_type(result_up_no_log, false);
				}
				else {
					return return_type(log2(result_up_no_log), true);
				}
			}
		}
		else {
			return return_type(up1+1, true);
		}
	}
	else { // (!upl1 && upl2)
		return autoint_add_upper_bound(up2, upl2, up1, upl1);
	}
}

constexpr std::pair<const int, const bool> autoint_mult_upper_bound(const int up1, const bool upl1, const int up2, const bool upl2) {
	using return_type = std::pair<const int, const bool>;
	if (upl1 && upl2) {
		return return_type(up1+up2, true);
	}
	else if (!upl1 && !upl2) {
		if (int64_t(up1)*int64_t(up2) <= std::numeric_limits<int32_t>::max()) { //fits int
			return return_type(up1*up2, false);
		}
		else {
			return return_type(log2(uint64_t(up1)*uint64_t(up2)), true);
		}
	}
	else if (upl1 && !upl2) {
		if (up1 < 63) { //up2 is at most 2^31-1 thus we can detect an overflow
			const uint64_t up1_no_log = (static_cast<uint64_t>(1) << up1)-1;
			const uint64_t result_up_no_log = up1_no_log * uint64_t(up2);
			if (result_up_no_log < up1_no_log) { //overflow
				return return_type(up1+log2(up2), true);
			}
			else {
				if (result_up_no_log <= std::numeric_limits<int32_t>::max()) {
					return return_type(result_up_no_log, false);
				}
				else {
					return return_type(log2(result_up_no_log), true);
				}
			}
		}
		else {
			return return_type(up1+log2(up2), true);
		}
	}
	else { // (!upl1 && upl2)
		return autoint_mult_upper_bound(up2, upl2, up1, upl1);
	}
}

///Upper bound is an absolute value
template<int T_UPPER_BOUND, bool T_UB_LOG = true, template<int> typename T_UNDERLYING_TYPE = AINT>
class AutoInt {
public:
	static constexpr int digits = (T_UB_LOG ? T_UPPER_BOUND : log2(T_UPPER_BOUND));
	using underlying_type = T_UNDERLYING_TYPE<digits>;
public:
	AutoInt() {}
	AutoInt(AutoInt && other) = default;
	AutoInt(AutoInt const & other) = default;
	
	AutoInt(underlying_type const & v) : m_v(v) {}
	AutoInt(underlying_type && v) : m_v(std::move(v)) {}
public:
	template<int T_O_UPPER_BOUND, bool T_O_UB_LOG>
	AutoInt & operator=(
		typename std::enable_if<
			(T_O_UB_LOG == T_UB_LOG && T_O_UPPER_BOUND <= T_UPPER_BOUND) ||
			(T_O_UB_LOG && !T_UB_LOG && T_O_UPPER_BOUND < log2(T_UPPER_BOUND)) ||
			(!T_O_UB_LOG && T_UB_LOG && log2(T_O_UPPER_BOUND) <= T_UPPER_BOUND),
			AutoInt<T_O_UPPER_BOUND, T_O_UB_LOG, T_UNDERLYING_TYPE>
		>::type const & other)
	{
		m_v = other.m_v;
		return *this;
	}
	AutoInt & operator=(AutoInt && other) {
		m_v = std::move(other.m_v);
		return *this;
	}
public:
	template<int T_O_UPPER_BOUND, bool T_O_UB_LOG>
	friend auto
	operator+(AutoInt<T_UPPER_BOUND, T_UB_LOG, T_UNDERLYING_TYPE> const & me, AutoInt<T_O_UPPER_BOUND, T_O_UB_LOG, T_UNDERLYING_TYPE> const & other) {
		constexpr std::pair<const int, const bool> result_traits = autoint_add_upper_bound(T_UPPER_BOUND, T_UB_LOG, T_O_UPPER_BOUND, T_O_UB_LOG);
		using return_type = AutoInt<result_traits.first, result_traits.second, T_UNDERLYING_TYPE>;
		using return_underlying_type = typename return_type::underlying_type;
		return return_type( ipt_static_cast<return_underlying_type>(me.value()) + ipt_static_cast<return_underlying_type>(other.value()) );
	}

	template<int T_O_UPPER_BOUND, bool T_O_UB_LOG>
	friend auto
	operator-(AutoInt<T_UPPER_BOUND, T_UB_LOG, T_UNDERLYING_TYPE> const & me, AutoInt<T_O_UPPER_BOUND, T_O_UB_LOG, T_UNDERLYING_TYPE> const & other) {
		constexpr std::pair<const int, const bool> result_traits = autoint_add_upper_bound(T_UPPER_BOUND, T_UB_LOG, T_O_UPPER_BOUND, T_O_UB_LOG);
		using return_type = AutoInt<result_traits.first, result_traits.second, T_UNDERLYING_TYPE>;
		using return_underlying_type = typename return_type::underlying_type;
		return return_type( ipt_static_cast<return_underlying_type>(me.value()) - ipt_static_cast<return_underlying_type>(other.value()) );
	}
	
	template<int T_O_UPPER_BOUND, bool T_O_UB_LOG>
	friend auto
	operator*(AutoInt<T_UPPER_BOUND, T_UB_LOG, T_UNDERLYING_TYPE> const & me, AutoInt<T_O_UPPER_BOUND, T_O_UB_LOG, T_UNDERLYING_TYPE> const & other) {
		constexpr std::pair<const int, const bool> result_traits = autoint_mult_upper_bound(T_UPPER_BOUND, T_UB_LOG, T_O_UPPER_BOUND, T_O_UB_LOG);
		using return_type = AutoInt<result_traits.first, result_traits.second, T_UNDERLYING_TYPE>;
		using return_underlying_type = typename return_type::underlying_type;
		return return_type( ipt_static_cast<return_underlying_type>(me.value()) * ipt_static_cast<return_underlying_type>(other.value()) );
	}
public:
	underlying_type const & value() const { return m_v; }
	underlying_type & value() { return m_v; }
private:
	
private:
	underlying_type m_v;
};

} // LIB_DTS2_NAMESPACE::detail::Kernel_sp

namespace LIB_RATSS_NAMESPACE {
	
template<int T_UPPER_BOUND, bool T_UB_LOG, template<int> typename T_UNDERLYING_TYPE>
struct Conversion< LIB_DTS2_NAMESPACE::detail::Kernel_sp::AutoInt<T_UPPER_BOUND, T_UB_LOG, T_UNDERLYING_TYPE> > {
	using type = LIB_DTS2_NAMESPACE::detail::Kernel_sp::AutoInt<T_UPPER_BOUND, T_UB_LOG, T_UNDERLYING_TYPE>;
	using underlying_type = typename type::underlying_type;
	static type moveFrom(mpq_class && v) {
		return type( Conversion<underlying_type>::moveFrom( std::move(v) ) );
	}
	static type moveFrom(mpq_class const & v) {
		return type( Conversion<underlying_type>::moveFrom(v) );
	}
	static mpq_class toMpq(const type & v) {
		return Conversion<underlying_type>::toMpq( v.value() );
	}
	static mpfr::mpreal toMpreal(const type & v, int precision) {
		return Conversion<underlying_type>::toMpreal( v.value(), precision);
	}
};

}

#endif
