#pragma once
#ifndef LIB_DTS2_KERNEL_SP_H
#define LIB_DTS2_KERNEL_SP_H

#include <libdts2/constants.h>
#include <libdts2/Kernel_sp/Point_sp.h>
#include <libdts2/Constrained_delaunay_triangulation_with_inexact_intersections_traits_s2.h>
#include <boost/multiprecision/cpp_int.hpp>

namespace LIB_DTS2_NAMESPACE {
namespace detail::Kernel_sp {

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
	static constexpr int bits = (requested_bits/64 + (requested_bits%64 == 0 ? 0 : 1))*64;
	using type = typename IntegerTypeFromBits<bits>::type;
};

template<typename T_BASE_TRAIT>
class Side_of_oriented_circle_s2:
	public T_BASE_TRAIT::Side_of_oriented_circle_2
{
public:
	using MyBaseTrait = T_BASE_TRAIT;
	using MyBaseClass = typename MyBaseTrait::Side_of_oriented_circle_2;
	using Is_auxiliary_point = typename MyBaseTrait::Is_auxiliary_point;
	using Point = Point_sp<T_BASE_TRAIT>;
	static constexpr int max_exponent = 30;
public:
	Side_of_oriented_circle_s2(MyBaseClass const & _base, Is_auxiliary_point const & _iap) :
	MyBaseClass(_base),
	m_iap(_iap)
	{}
public:
	CGAL::Sign operator()(Point const & p, Point const & q, Point const & r, Point const & t) {
		if (p.pos == q.pos() && p.pos() == r.pos() && p.pos() == t.pos() &&
			!m_iap(p) && !m_iap(q) && !m_iap(r) && !m_iap(t) &&
			p.exponent() < max_exponent && q.exponent() < max_exponent && r.exponent() < max_exponent && t.exponent() < max_exponent
		)
		{
			//we have to normalize the denominators
			int64_t max_den = std::max(std::max(p.denominator(), q.denominator()), std::max(r.denominator(), t.denominator()));
			return calc<30>(
				p.numerator0()*p.denominator()/max_den, p.numerator1()*p.denominator()/max_den,
				q.numerator0()*q.denominator()/max_den, q.numerator1()*q.denominator()/max_den,
				r.numerator0()*r.denominator()/max_den, r.numerator1()*r.denominator()/max_den,
				t.numerator0()*t.denominator()/max_den, t.numerator1()*t.denominator()/max_den
			);
		}
		else { //use base predicate
			return MyBaseClass::operator()(p, q, r, t);
		}
	}
public:
	//T_STAGE0: n Bits
	//T_STAGE1: n+1 Bits
	//T_STAGE2: 2n+2 Bits
	//T_STAGE3: 2n+3 Bits
	//T_STAGE4: 4n+6 Bits
	template<
		int T_START_BITS,
		typename T_STAGE0 = typename AlignedIntegerTypeFromBits<T_START_BITS>::type,
		typename T_STAGE1 = typename AlignedIntegerTypeFromBits<T_START_BITS+1>::type,
		typename T_STAGE2 = typename AlignedIntegerTypeFromBits<2*T_START_BITS+2>::type,
		typename T_STAGE3 = typename AlignedIntegerTypeFromBits<2*T_START_BITS+3>::type,
		typename T_STAGE4 = typename AlignedIntegerTypeFromBits<4*T_START_BITS+6>::type
	>
	inline CGAL::Sign calc(
					const T_STAGE0 & p_num0, const T_STAGE0 & p_num1,
					const T_STAGE0 & q_num0, const T_STAGE0 & q_num1,
					const T_STAGE0 & r_num0, const T_STAGE0 & r_num1,
					const T_STAGE0 & t_num0, const T_STAGE0 & t_num1
	)
	{
		//points have up to T_STAGE0=n Bits
		T_STAGE1 px(p_num0), py(p_num1);
		T_STAGE1 qx(q_num0), qy(q_num1);
		T_STAGE1 rx(r_num0), ry(r_num1);
		T_STAGE1 tx(t_num0), ty(t_num1);
		//now up to n Bits
		
		T_STAGE2 qpx = qx-px;
		T_STAGE2 qpy = qy-py;
		T_STAGE2 rpx = rx-px;
		T_STAGE2 rpy = ry-py;
		T_STAGE2 tpx = tx-px;
		T_STAGE2 tpy = ty-py;
		//now up to n+1 Bits
		
		//intermediates are 2*(n+1) = 2n+2 Bits and another Bit for the addition
		T_STAGE4 a00 = T_STAGE3(qpx*tpy) - T_STAGE3(qpy*tpx);
		T_STAGE4 a01 = T_STAGE3(tpx*T_STAGE2(tx-qx)) + T_STAGE3(tpy*T_STAGE2(ty-qy));
		T_STAGE4 a10 = T_STAGE3(qpx*rpy) - T_STAGE3(qpy*rpx);
		T_STAGE4 a11 = T_STAGE3(rpx*T_STAGE2(rx-qx)) + T_STAGE3(rpy*T_STAGE2(ry-qy));
		//now up to 2n+3
		
		//determinant is a00*a11 - a10*a01
		T_STAGE4 a0011 = a00*a11;
		T_STAGE4 a1001 = a10*a01;
		//now up to 2*(2n+3) = 4n+6 Bits
		
		a1001 = -a1001;
		
		//TODO: get rid of the following branch cascade
		//Simply do a0011 - a1001 and check if a signed overflow occured
		
		if (a0011 < 0 && a1001 < 0) {
			return CGAL::NEGATIVE;
		}
		else if (a0011 > 0 && a1001 > 0) {
			return CGAL::POSITIVE;
		}
		else if (a0011 == 0) {
			return a1001 < 0 ? CGAL::NEGATIVE : (a1001 > 0 ? CGAL::POSITIVE : CGAL::ZERO);
		}
		else if (a1001 == 0) {
			return a0011 < 0 ? CGAL::NEGATIVE : (a0011 > 0 ? CGAL::POSITIVE : CGAL::ZERO);
		}
		else if (a0011 < 0) { // && a1001 > 0
			a0011 = -a0011;
			if (a0011 > a1001) {
				return CGAL::NEGATIVE;
			}
			else if (a0011 == a1001) {
				return CGAL::NEGATIVE;
			}
			else {
				return CGAL::POSITIVE;
			}
		}
		else { // a0011 > 0 && a1001 < 0
			a1001 = -a1001;
			if (a1001 > a0011) {
				return CGAL::NEGATIVE;
			}
			else if (a0011 == a1001) {
				return CGAL::ZERO;
			}
			else { //a1001 < a0011
				return CGAL::POSITIVE;
			}
		}
	}
private:
	Is_auxiliary_point m_iap;
};

template<typename T_BASE_TRAIT>
class Orientation_s2: public T_BASE_TRAIT::Orientation_2 {
public:
	using MyBaseTrait = T_BASE_TRAIT;
	using MyBaseClass = typename MyBaseTrait::Orientation_2;
	using Is_auxiliary_point = typename MyBaseTrait::Is_auxiliary_point;
	using Point = Point_sp<T_BASE_TRAIT>;
	static constexpr int max_exponent = 30;
public:
	CGAL::Sign operator()(Point const & p, Point const & q, Point const & t) {
		if (p.pos == q.pos() && p.pos() == t.pos()) {
			assert(p.exponent() == q.exponent() && p.exponent() == t.exponent());
			//calculate a point reflected third point,
			//scale all points by the common denominator
			//and use the Oriented_side_of_circle predicate
			int128 r_num0 = -int64(p.numerator0()) * int64(p.denominator());
			int128 r_num1 = -int64(p.numerator1()) * int64(p.denominator());
			int128 r_den = int64(p.numerator0())*int64(p.numerator0()) - int64(p.numerator1())*int64(p.numerator1());
			
// 			int128 com_den = r_den * p.denominator();

			//now go to a common denominator and forget about the denominator
			int128 p_num0 = p.numerator0() * r_den;
			int128 p_num1 = p.numerator1() * r_den;
			int128 q_num0 = q.numerator0() * r_den;
			int128 q_num1 = q.numerator1() * r_den;
			r_num0 *= p.denominator();
			r_num1 *= p.denominator();
			int128 t_num0 = t.numerator0() * r_den;
			int128 t_num1 = t.numerator1() * r_den;
			//all points now have less than 92 Bits
			
			return m_soc.calc<92>(
					p_num0, p_num1,
					q_num0, q_num1,
					r_num0, r_num1,
					t_num0, t_num1
			);
		}
		else { //use base kernel
			return MyBaseClass::operator()(p, q, t);
		}
	}
private:
	Side_of_oriented_circle_s2<BaseTrait> m_soc;
};


} //end namespace detail::Kernel_sp 


template<
	typename T_LINEAR_KERNEL,
	typename T_AUX_POINT_GENERATOR=detail::EpsBasedAuxPoints<T_LINEAR_KERNEL>
>
class Kernel_sp:
	public Constrained_delaunay_triangulation_with_inexact_intersections_base_traits_s2<
		T_LINEAR_KERNEL,
		T_AUX_POINT_GENERATOR,
		Point_sp<T_LINEAR_KERNEL>
		>
{
public:
	static constexpr int snap_bits = 31;
public:
	using Self = Kernel_sp<T_LINEAR_KERNEL, T_AUX_POINT_GENERATOR>;
	using MyBaseTrait = Constrained_delaunay_triangulation_with_inexact_intersections_base_traits_s2<
		T_LINEAR_KERNEL,
		T_AUX_POINT_GENERATOR,
		Point_sp<T_LINEAR_KERNEL>
	>;
	using AuxiliaryPointsGenerator = typename MyBaseTrait::AuxiliaryPointsGenerator;
public:
	using Side_of_oriented_circle_2 = detail::Kernel_sp::Side_of_oriented_circle_s2<MyBaseTrait>;
	using Orientation_2 = detail::Kernel_sp::Orientation_s2<MyBaseTrait>;
public:
	Kernel_sp();
	Kernel_sp(AuxiliaryPointsGenerator const & _apg, int _significands) :
	MyBaseTrait(_apg, _significands)
	{}
	Kernel_sp(const Kernel_sp & other) :
	MyBaseTrait(other)
	{}
	~Kernel_sp() override {}
public:
	Side_of_oriented_circle_2 side_of_oriented_circle_2_object () const;
	Orientation_2 orientation_2_object () const;
};

}//end namespace

#endif
