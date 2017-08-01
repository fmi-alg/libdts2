#ifndef LIB_DTS2_KERNEL_SP_H
#define LIB_DTS2_KERNEL_SP_H
#pragma once

#include <libdts2/constants.h>
#include <libdts2/Point_sp.h>
#include <libdts2/Constrained_delaunay_triangulation_with_inexact_intersections_traits_s2.h>
#include <boost/multiprecision/cpp_int.hpp>

namespace LIB_DTS2_NAMESPACE {
namespace internal {
namespace Kernel_sp {

template<typename T_BASE_TRAIT>
class Oriented_side_of_circle_s2{
public:
	using BaseTrait = T_BASE_TRAIT;
	using Point = Point_sp<T_BASE_TRAIT>;
	using Point_2 = Point_sp<T_BASE_TRAIT>;
	using Point_3 = Point_sp<T_BASE_TRAIT>;
	using int128 = __int128_t;
	using int64 = int64_t;
public:
	CGAL::Sign operator()(Point const & p, Point const & q, Point const & r, Point const & t) {
		if (p.pos == q.pos() && p.pos() == r.pos() && p.pos() == t.pos() &&
			!p.is_auxiliary() && !q.is_auxiliary() && !r.is_auxiliary() && !t.is_auxiliary()
		)
		{
		
			//points have up to 30 Bits
			int64_t px(p.num0()), py(p.num1());
			int64_t qx(q.num0()), qy(p.num1());
			int64_t rx(r.num0()), ry(p.num1());
			int64_t tx(t.num0()), ty(t.num1());
			//now up to 30 Bits
			
			int64_t qpx = qx-px;
			int64_t qpy = qy-py;
			int64_t rpx = rx-px;
			int64_t rpy = ry-py;
			int64_t tpx = tx-px;
			int64_t tpy = ty-py;
			//now up to 30+1 Bits
			
			int128 a00 = int128(int64(qpx)*int64(tpy)) - int128(int64(qpy)*int64(tpx));
			int128 a01 = int128(int64(tpx)*int64(tx-qx)) + int128(int64(tpy)*int64(ty-qy));
			int128 a10 = int128(int64(qpx)*int64(rpy)) - int128(int64(qpy)*int64(rpx));
			int128 a11 = int128(int64(rpx)*int64(rx-qx)) + int128(int64(rpy)*int64(ry-qy));
			//now up to (30+1)*2+1 = 63 Bits
			
			//determinant is a00*a11 - a10*a01
			int128 a0011 = a00*a11;
			int128 a1001 = a10*a01;
			//now up to ((30+1)*2+1)*2 = 126 Bits
			
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
				return sign(a1001);
			}
			else if (a1001 == 0) {
				return sign(a0011);
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
		else { //use base predicate
			return (*this)(p.point3(), q.point3(), r.point3(), t.point3());
		}
	}
public:
	//T_STAGE0: n Bits
	//T_STAGE1: n+1 Bits
	//T_STAGE2: 2n+2 Bits
	//T_STAGE3: 2n+3 Bits
	//T_STAGE4: 4n+6 Bits
	template<typename T_STAGE0, typename T_STAGE1, typename T_STAGE2, typename T_STAGE3, typename T_STAGE4>
	CGAL::Sign calc(const T_STAGE0 & p_num0, const T_STAGE0 & p_num1,
					const T_STAGE0 & q_num0, const T_STAGE0 & q_num1,
					const T_STAGE0 & r_num0, const T_STAGE0 & r_num1,
					const T_STAGE0 & t_num0, const T_STAGE0 & t_num1
	)
	{
		//points have up to T_STAGE0 Bits
		T_STAGE1 px(p_num0), py(p_num1);
		T_STAGE1 qx(q_num0), qy(q_num1);
		T_STAGE1 rx(r_num0), ry(r_num1);
		T_STAGE1 tx(t_num0), ty(t.num1);
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
			return sign(a1001);
		}
		else if (a1001 == 0) {
			return sign(a0011);
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
	CGAL::Sign operator()(const Point_3 & p, const Point_3 & q, const Point_3 & r, const Point_3 & t) {
		return m_bp(p, q, r, t);
	}
private:
	CGAL::Sign sign(int128 v) const {
		if (v < 0) {
			return CGAL::NEGATIVE;
		}
		else if (v == 0) {
			return CGAL::ZERO;
		}
		else {
			return CGAL::POSITIVE;
		}
	}
public:
	BaseTrait::Oriented_side_of_circle_2 m_bp;
};

template<typename T_BASE_TRAIT>
class Orientation_s2 {
public:
	using BaseTrait = T_BASE_TRAIT;
	using LinearKernel = typename BaseTrait::LinearKernel;
	using Point = Point_sp;
	using Point_2 = Point_sp;
	using Point_3 = Point_sp;
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
			
			//T_STAGE0: n Bits; 92
			//T_STAGE1: n+1 Bits; 93
			//T_STAGE2: 2n+2 Bits; 186
			//T_STAGE3: 2n+3 Bits; 187
			//T_STAGE4: 4n+6 Bits; 374
			return m_soc.calc<int128, int128, boost_int<186>, boost_int<187>, boost_int<374>>(
					p_num0, p_num1,
					q_num0, q_num1,
					r_num0, r_num1,
					t_num0, t_num1
			);
		}
		else { //use base kernel
			(*this)(p.point3(), q.point3(), t.point3());
		}
	}
private:
	template<int T_V>
	struct Align64 {
		static constexpr const int result = (T_V/64 + (T_V%64 == 0 ? 0 : 1))*64; 
	};
	template<int T_BITS>
	using boost_int = boost::multiprecision::number<boost::multiprecision::cpp_int_backend<
		Align64<T_BITS>::result, Align64<T_BITS>::result, signed_magnitude, unchecked, void> >;
	using int64 = int64_t;
	using int128 = __int128;
private:
	CGAL::Sign operator()(const Point_3 & p, const Point_3 & q, const Point_3 & r) {
		return m_ot2(p, q, r);
	}
	
private:
	BaseTrait::Orientation_2 m_ot2;
	Oriented_side_of_circle_s2<BaseTrait> m_soc;
};


}}


template<typename T_BASE_TRAIT>
class Kernel_sp: 

}//end namespace

#endif