#pragma once
#ifndef LIB_DTS2_KERNEL_SP_H
#define LIB_DTS2_KERNEL_SP_H

#include <libdts2/constants.h>
#include <libdts2/Kernel_sp/Point_sp.h>
#include <libdts2/Constrained_delaunay_triangulation_with_inexact_intersections_traits_s2.h>
#include <libratss/Conversion.h>

namespace LIB_DTS2_NAMESPACE {
namespace detail::Kernel_sp {

template<typename T_BASE_TRAIT>
class Side_of_oriented_circle_s2:
	protected T_BASE_TRAIT::Side_of_oriented_circle_2
{
public:
	using MyBaseTrait = T_BASE_TRAIT;
	using MyBaseClass = typename MyBaseTrait::Side_of_oriented_circle_2;
	using Is_auxiliary_point = typename MyBaseTrait::Is_auxiliary_point;
	using Point = Point_sp<typename MyBaseTrait::LinearKernel>;
	static constexpr int max_exponent = 31;
public:
	Side_of_oriented_circle_s2(MyBaseClass const & _base, Is_auxiliary_point const & _iap) :
	MyBaseClass(_base),
	m_iap(_iap)
	{}
	MyBaseClass const & base() const { return *this; }
	MyBaseClass & base() { return *this; }
public:
	CGAL::Sign operator()(Point const & p, Point const & q, Point const & r, Point const & t) const {
		if (p.pos() == q.pos() && p.pos() == r.pos() && p.pos() == t.pos() &&
			p.exponent() < max_exponent && q.exponent() < max_exponent && r.exponent() < max_exponent && t.exponent() < max_exponent &&
			!m_iap(p) && !m_iap(q) && !m_iap(r) && !m_iap(t)
		)
		{
			//we have to normalize the denominators
			auto max_den = std::max(std::max(p.denominator(), q.denominator()), std::max(r.denominator(), t.denominator()));
			int result = calc<max_exponent>(
				p.numerator0()*(max_den/p.denominator()), p.numerator1()*(max_den/p.denominator()),
				q.numerator0()*(max_den/q.denominator()), q.numerator1()*(max_den/q.denominator()),
				r.numerator0()*(max_den/r.denominator()), r.numerator1()*(max_den/r.denominator()),
				t.numerator0()*(max_den/t.denominator()), t.numerator1()*(max_den/t.denominator())
			);
			switch (std::abs<int>(p.pos())) {
			case 1:
			case 3:
				result *= (p.pos() < 0 ? -1 : 1);
				break;
			case 2: //This due to the fact that we store num0/num1 such that handedness of the coordinate system changes if we project through the y-axis
				result *=  (p.pos() < 0 ? 1 : -1);
				break;
			default:
				throw std::runtime_error("Invalid position on sphere: " + std::to_string(p.pos()));
			};
			#ifndef NDEBUG
			{
				using SOC = typename MyBaseTrait::LinearKernel::Side_of_oriented_circle_2;
				using FTLK = typename MyBaseTrait::LinearKernel::FT;
				using Point_2LK = typename MyBaseTrait::LinearKernel::Point_2;
				SOC soc;
				
				auto conv = [](mpq_class const & x) -> FTLK { return ratss::convert<FTLK>(x); };
				
				#define ASSIGN_POINT(__NAME) Point_2LK __NAME ## lk(conv(__NAME.numerator0())*conv(max_den/__NAME.denominator()), conv(__NAME.numerator1())*conv(max_den/__NAME.denominator()));
				ASSIGN_POINT(p)
				ASSIGN_POINT(q)
				ASSIGN_POINT(r)
				ASSIGN_POINT(t)
				#undef ASSIGN_POINT
				
				auto should = MyBaseClass::operator()(p, q, r, t);
				int should_2d;
				switch (std::abs<int>(p.pos())) {
				case 1:
				case 3:
					should_2d = ((p.pos() < 0 ? -1 : 1) * soc(plk, qlk, rlk, tlk));
					break;
				case 2: 
					should_2d =  ((p.pos() < 0 ? 1 : -1) * soc(plk, qlk, rlk, tlk));
					break;
				default:
					throw std::runtime_error("Invalid position on sphere: " + std::to_string(p.pos()));
				};
				assert(should == should_2d);
				assert(should_2d == result);
				assert(result == should);
			}
			#endif
			return CGAL::Sign(result);
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
	) const
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
public:
	Is_auxiliary_point const & iap() const { return m_iap; }
private:
	Is_auxiliary_point m_iap;
};

template<typename T_BASE_TRAIT>
class Orientation_s2: public T_BASE_TRAIT::Orientation_2 {
public:
	using MyBaseTrait = T_BASE_TRAIT;
	using MyBaseClass = typename MyBaseTrait::Orientation_2;
	using Is_auxiliary_point = typename MyBaseTrait::Is_auxiliary_point;
	using Point = Point_sp<typename MyBaseTrait::LinearKernel>;
	static constexpr int max_exponent = 30;
public:
	Orientation_s2(MyBaseClass const & _base, Side_of_oriented_circle_s2<MyBaseTrait> const & _soc) :
	MyBaseClass(_base),
	m_soc(_soc)
	{}
public:
	CGAL::Sign operator()(Point const & p, Point const & q, Point const & t) const {
		return MyBaseClass::operator()(p, q, t);
		if (p.pos() == q.pos() && p.pos() == t.pos() &&
			!m_soc.iap()(p) && !m_soc.iap()(q) && !m_soc.iap()(t) &&
			p.exponent() < max_exponent && q.exponent() < max_exponent && t.exponent() < max_exponent
		)
		{
			//we have to normalize the denominators
			auto max_den = std::max(std::max(p.denominator(), q.denominator()), std::max(t.denominator(), t.denominator()));
			auto result = (CGAL::Sign) -calc<30>(
				p.numerator0()*(max_den/p.denominator()), p.numerator1()*(max_den/p.denominator()),
				q.numerator0()*(max_den/q.denominator()), q.numerator1()*(max_den/q.denominator()),
				t.numerator0()*(max_den/t.denominator()), t.numerator1()*(max_den/t.denominator())
			);
			#ifndef NDEBUG
			{
				auto should = MyBaseClass::operator()(p, q, t);
				if (should != result) {
					should = MyBaseClass::operator()(p, q, t);
					result = (CGAL::Sign) -calc<30>(
						p.numerator0()*(max_den/p.denominator()), p.numerator1()*(max_den/p.denominator()),
						q.numerator0()*(max_den/q.denominator()), q.numerator1()*(max_den/q.denominator()),
						t.numerator0()*(max_den/t.denominator()), t.numerator1()*(max_den/t.denominator())
					);
				}
				assert(result == should);
			}
			#endif
			return result;
		}
		else { //use base predicate
			return MyBaseClass::operator()(p, q, t);
		}
	}
	template<
		int T_START_BITS,
		typename T_STAGE0 = typename AlignedIntegerTypeFromBits<T_START_BITS>::type,
		typename T_STAGE1 = typename AlignedIntegerTypeFromBits<T_START_BITS+1>::type
	>
	CGAL::Sign calc(
					const T_STAGE0 & p_num0, const T_STAGE0 & p_num1,
					const T_STAGE0 & q_num0, const T_STAGE0 & q_num1,
					const T_STAGE0 & t_num0, const T_STAGE0 & t_num1
	) const {
		//calculate a point reflected third point
		//and use the Oriented_side_of_circle predicate
		T_STAGE1 r_num0 = -T_STAGE1(p_num0);
		T_STAGE1 r_num1 = -T_STAGE1(p_num1);
		return m_soc.template calc<T_START_BITS+1, T_STAGE1> (
				p_num0, p_num1,
				q_num0, q_num1,
				r_num0, r_num1,
				t_num0, t_num1
		);
	}
private:
	Side_of_oriented_circle_s2<MyBaseTrait> m_soc;
};

#define LESS_VAR_THREE(__VAR) template<typename T_BASE_TRAIT> \
class Less_ ## __VAR ## _3: public T_BASE_TRAIT::Less_ ## __VAR ## _3 { \
public: \
	using MyBaseTrait = T_BASE_TRAIT; \
	using MyParent = typename MyBaseTrait::Less_ ## __VAR ## _3; \
	using Point = Point_sp<typename MyBaseTrait::LinearKernel>; \
public: \
	Less_ ## __VAR ## _3(); \
	Less_ ## __VAR ## _3(MyParent const & other): MyParent(other) {} \
public: \
	inline bool operator()(Point const & a, Point const & b) const { \
		return a.template __VAR <mpq_class>() < b.template __VAR <mpq_class>(); \
	} \
};

LESS_VAR_THREE(x)
LESS_VAR_THREE(y)
LESS_VAR_THREE(z)

#undef LESS_VAR_TREE

template<typename T_BASE_TRAIT>
class Construct_point_sp: public T_BASE_TRAIT::Construct_point_2 {
public:
	using MyBaseTrait = T_BASE_TRAIT;
	using MyBaseClass = typename MyBaseTrait::Construct_point_2;
	using FT = typename MyBaseTrait::FT;
	using Point = Point_sp<typename MyBaseTrait::LinearKernel>;
public:
	Construct_point_sp() {}
	Construct_point_sp(MyBaseClass const & base) : MyBaseClass(base) {}
public:
	template<typename ...T>
	Point operator()(T && ... t) const { return Point(std::forward<T>(t)...); }
};

} //end namespace detail::Kernel_sp 


template<
	typename T_LINEAR_KERNEL,
	typename T_AUX_POINT_GENERATOR=detail::EpsBasedAuxPoints<T_LINEAR_KERNEL>
>
class Kernel_sp_base:
	public Constrained_delaunay_triangulation_with_inexact_intersections_base_traits_s2<
		T_LINEAR_KERNEL,
		T_AUX_POINT_GENERATOR,
		Point_sp<T_LINEAR_KERNEL>
		>
{
public:
	static constexpr int snap_bits = 31;
public:
	using Self = Kernel_sp_base<T_LINEAR_KERNEL, T_AUX_POINT_GENERATOR>;
	using MyBaseTrait = Constrained_delaunay_triangulation_with_inexact_intersections_base_traits_s2<
		T_LINEAR_KERNEL,
		T_AUX_POINT_GENERATOR,
		Point_sp<T_LINEAR_KERNEL>
	>;
	using AuxiliaryPointsGenerator = typename MyBaseTrait::AuxiliaryPointsGenerator;
public:
	using Less_x_3 = detail::Kernel_sp::Less_x_3<MyBaseTrait>;
	using Less_y_3 = detail::Kernel_sp::Less_y_3<MyBaseTrait>;
	using Less_z_3 = detail::Kernel_sp::Less_z_3<MyBaseTrait>;
	
	using Less_x_2 = Less_x_3;
	using Less_y_2 = Less_y_3;
	
	using Construct_point_2 = detail::Kernel_sp::Construct_point_sp<MyBaseTrait>;
public:
	using Side_of_oriented_circle_2 = detail::Kernel_sp::Side_of_oriented_circle_s2<MyBaseTrait>;
	using Orientation_2 = detail::Kernel_sp::Orientation_s2<MyBaseTrait>;
public:
	Kernel_sp_base() {}
	Kernel_sp_base(AuxiliaryPointsGenerator const & _apg, int _significands, int _intersectSignificands = -1) :
	MyBaseTrait(_apg, _significands, _intersectSignificands)
	{}
	Kernel_sp_base(const Kernel_sp_base & other) :
	MyBaseTrait(other)
	{}
	~Kernel_sp_base() {}
public:
	Side_of_oriented_circle_2 side_of_oriented_circle_2_object () const {
		return Side_of_oriented_circle_2( MyBaseTrait::side_of_oriented_circle_2_object(), MyBaseTrait::is_auxiliary_point_object() );
	}
	Orientation_2 orientation_2_object() const {
		return Orientation_2( MyBaseTrait::orientation_2_object(), Self::side_of_oriented_circle_2_object() );
	}
	
	#define LESS_VAR_OBJECT(__VAR) Less_ ## __VAR less_ ## __VAR ## _object() const { return Less_ ## __VAR( MyBaseTrait::less_ ## __VAR ## _object() ); }
	LESS_VAR_OBJECT(x_3)
	LESS_VAR_OBJECT(y_3)
	LESS_VAR_OBJECT(z_3)
	LESS_VAR_OBJECT(x_2)
	LESS_VAR_OBJECT(y_2)
	#undef LESS_VAR_OBJECT
	
	Construct_point_2 construct_point_2_object() const { return Construct_point_2( MyBaseTrait::construct_point_2_object() ); }
	
};

using Kernel_sp = Kernel_sp_base<CGAL::Exact_predicates_exact_constructions_kernel>;
using Kernel_sp_64 = Kernel_sp_base<CGAL::Filtered_simple_cartesian_extended_integer_kernel>;

}//end namespace

#endif
