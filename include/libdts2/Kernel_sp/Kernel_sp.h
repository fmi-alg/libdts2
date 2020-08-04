#pragma once
#ifndef LIB_DTS2_KERNEL_SP_H
#define LIB_DTS2_KERNEL_SP_H

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
			return calc<30>();
		}
		else { //use base predicate
			return (*this)(p.point3(), q.point3(), r.point3(), t.point3());
		}
	}
public:
	template<int T_START_BITS>
	CGAL::Sign calc(const T_STAGE0 & p_num0, const T_STAGE0 & p_num1,
					const T_STAGE0 & q_num0, const T_STAGE0 & q_num1,
					const T_STAGE0 & r_num0, const T_STAGE0 & r_num1,
					const T_STAGE0 & t_num0, const T_STAGE0 & t_num1
	)
	{
		calc<
			boost_int<T_START_BITS>,
			boost_int<T_START_BITS>,
			boost_int<T_START_BITS>,
			boost_int<T_START_BITS>,
			boost_int<T_START_BITS>
		>();
	}
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
	CGAL::Sign operator()(const Point_3 & p, const Point_3 & q, const Point_3 & r) {
		return m_ot2(p, q, r);
	}
	
private:
	BaseTrait::Orientation_2 m_ot2;
	Oriented_side_of_circle_s2<BaseTrait> m_soc;
};


}}


template<typename T_BASE_TRAIT>
class Kernel_sp {
public:
	static constexpr int snap_bits = 30;
public:
	using Self = Kernel_sp<T_BASE_TRAIT>;
	using BaseTraits = T_BASE_TRAIT;
	using LinearKernel = typename BaseTraits::LinearKernel;
	using FT = typename BaseTraits::FT;
	using Projector = typename BaseTraits::Projector;
public:
	using Segment_2 
public:
	using dummy = int;
	using Side_of_oriented_circle_2 = dummy;
	using Orientation_2 = dummy;
	using Collinear_are_ordered_along_line_2 = dummy;
public:
	using Compare_distance_2 = dummy;
	using Compare_x_2 = dummy;
	using Compare_y_2 = dummy;
	using Less_x_2 = dummy;
	using Less_y_2 = dummy;
	using Less_x_3 = dummy;
	using Less_y_3 = dummy;
	using Less_z_3 = dummy;
public:
	using Construct_segment_2 = dummy;
	using Construct_line_2 = dummy;
	using Intersect_2 = dummy;
	using Compute_squared_distance_2 = dummy;
	using Project_on_sphere = dummy;
public: //accessor functions
	const FT & epsilon() const;
public:
	Kernel_sp();
	Kernel_sp(const FT & _epsilon) :
	m_traits(_epsilon, snap_bits)
	{}
	Kernel_sp(const Kernel_sp & other) :
	m_traits(other.m_traits)
	{}
	~Kernel_sp() {}
	const BaseTraits & baseTraits() const { return m_traits; }
public:
	int significands() const { return baseTraits().significands(); }
	int intersectSignificands() const { return significands(); }
	const FT & epsZ() const { return baseTraits().epsZ(); }
	const Projector & projector() const { return baseTraits().projector(); }
public: //objects
	Side_of_oriented_circle_2 side_of_oriented_circle_2_object () const;
	Orientation_2 orientation_2_object () const;
	Collinear_are_ordered_along_line_2 collinear_are_ordered_along_line_2_object() const;
	
	Compare_distance_2 compare_distance_2_object () const;
	Compare_x_2 compare_x_2_object() const ;
	Compare_y_2 compare_y_2_object() const;
	
	Less_x_2 less_x_2_object() const;
	Less_y_2 less_y_2_object() const;
	Less_x_3 less_x_3_object() const;
	Less_y_3 less_y_3_object() const;
	Less_z_3 less_z_3_object() const;

	Construct_segment_2 construct_segment_2_object() const;
	Intersect_2 intersect_2_object() const;
	Construct_line_2 construct_line_2_object() const;

	Compute_squared_distance_2 compute_squared_distance_2_object() const;
	Project_on_sphere project_on_sphere_object() const;
	Project_on_sphere project_on_sphere_object(int _significands) const;
private:
	BaseTraits m_traits;
};

}//end namespace

#endif
