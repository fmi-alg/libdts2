#ifndef LIB_DTS2_KERNEL_SP_H
#define LIB_DTS2_KERNEL_SP_H
#pragma once

#include <libdts2/constants.h>
#include <libdts2/Point_sp.h>

namespace LIB_DTS2_NAMESPACE {
namespace internal {
namespace Kernel_sp {

template<typename T_BASE_TRAIT>
class Orientation_s2 {
public:
	using BaseTrait = T_BASE_TRAIT;
	using Point = Point_sp;
	using Point_2 = Point_sp;
	using Point_3 = Point_sp;
public:
	CGAL::Sign operator()(Point const & p, Point const & q, Point const & r) {
		if (p.pos == q.pos() && p.pos() == r.pos()) {
			
		}
		else { //use base kernel
			
		}
	}
private:
	BaseTrait::Orientation_3 m_ot3;
};

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
			
			int128 det = a0011 - a1001;
			
			if (det < 0) {
				return CGAL::NEGATIVE;
			}
			else if (det == 0) {
				return CGAL::ZERO;
			}
			else {
				return CGAL::POSITIVE;
			}
		}
		else { //use base predicate
			return m_bp(p.point3(), q.point3(), r.point3(), t.point3());
		}
	}
public:
	BaseTrait::Oriented_side_of_circle_2 m_bp;
};

}}


class Kernel_sp {



};

}//end namespace

#endif