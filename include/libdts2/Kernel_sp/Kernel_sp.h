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
	static constexpr int max_exponent = 32;
	//Signed integer which can represent an unsigned number with up to N bits
	template<int N>
	using AINT = typename AlignedIntegerTypeFromBits<N+1>::type; //Plus one for the sign
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
	template<int N>
	inline CGAL::Sign calc(
					const AINT<N> & p_num0, const AINT<N> & p_num1,
					const AINT<N> & q_num0, const AINT<N> & q_num1,
					const AINT<N> & r_num0, const AINT<N> & r_num1,
					const AINT<N> & t_num0, const AINT<N> & t_num1
	) const
	{
		
		using T_STAGE1 = AINT<N+1>;
		using T_STAGE2 = AINT<2*N+2>;
		using T_STAGE3 = AINT<2*N+3>;
		using T_STAGE4 = AINT<4*N+6>;
		
		//points have up to AINT<N>=T_STAGE0=n Bits
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
	using Is_in_auxiliary_triangle = typename MyBaseTrait::Is_in_auxiliary_triangle;
	using Point = Point_sp<typename MyBaseTrait::LinearKernel>;
	static constexpr int max_exponent = 32;
	template<int N>
	using AINT = typename AlignedIntegerTypeFromBits<N+1>::type; //Plus 1 for the sign
private:
	using Numerators_3 = typename Point::Numerators_3;
public:
	Orientation_s2(MyBaseClass const & _base, Is_in_auxiliary_triangle const & _iat) :
	MyBaseClass(_base),
	m_iat(_iat)
	{}
public:
	CGAL::Sign operator()(Point const & p, Point const & q, Point const & t) const {
		auto max_exp = std::max({p.exponent(), q.exponent(), t.exponent()});
		if (max_exp < max_exponent) {
			Numerators_3 p3 = p.numerators3();
			Numerators_3 q3 = q.numerators3();
			Numerators_3 origin(LIB_DTS2_ORIGIN_X, LIB_DTS2_ORIGIN_Y, LIB_DTS2_ORIGIN_Z);
			Numerators_3 t3 = t.numerators3();
			
			CGAL::Sign oriented_side = ot3(p3, q3, origin, t3, max_exp);
			
			if (oriented_side != CGAL::Sign::COLLINEAR && m_iat(p) && m_iat(q)) {
				//the point always has to be on the opposite side of the infinite vertex
				oriented_side = - ot3(p3, q3, origin, Numerators_3(0, 0, 1), max_exp);
			}
			assert(oriented_side == MyBaseClass::operator()(p, q, t));
			return oriented_side;
		}
		else { //use base predicate
			return MyBaseClass::operator()(p, q, t);
		}
	}
private:
	CGAL::Sign ot3(Numerators_3 const & p, Numerators_3 const & q, Numerators_3 const & r, Numerators_3 const & s, uint32_t /*max_exp*/) const {
		return ot3<max_exponent>(p, q, r, s);
	}
	template<int N>
	CGAL::Sign ot3(Numerators_3 const & p, Numerators_3 const & q, Numerators_3 const & r, Numerators_3 const & s) const {
		return ot3<2*N+1>(p.x, p.y, p.z,
				   q.x, q.y, q.z,
				   r.x, r.y, r.z,
				   s.x, s.y, s.z
			   );
	}
	template<int N>
	CGAL::Sign ot3(
					const AINT<N+1> & px, const AINT<N+1> & py, const AINT<N+1> & pz,
					const AINT<N+1> & qx, const AINT<N+1> & qy, const AINT<N+1> & qz,
					const AINT<N+1> & rx, const AINT<N+1> & ry, const AINT<N+1> & rz,
					const AINT<N+1> & sx, const AINT<N+1> & sy, const AINT<N+1> & sz
	) const {
		return sign_of_determinant<N+1>(qx-px,rx-px,sx-px,
                                    qy-py,ry-py,sy-py,
                                    qz-pz,rz-pz,sz-pz);
	}
	
	template<int N>
	CGAL::Sign sign_of_determinant(const AINT<2*N> & a00,  const AINT<2*N>& a01,  const AINT<2*N>& a02,
									const AINT<2*N> & a10,  const AINT<2*N>& a11,  const AINT<2*N>& a12,
									const AINT<2*N> & a20,  const AINT<2*N>& a21,  const AINT<2*N>& a22) const
	{
		using T_STAGE1 = AINT<2*N+1>;
		using T_STAGE2 = AINT<3*N+1>;
		using T_STAGE3 = AINT<3*N+3>;
		//aij have 2*n bits
		//based on CGAL's sign_of_determinant
				// First compute the det2x2
		T_STAGE2 m01 = T_STAGE1(a00*a11) - T_STAGE1(a10*a01);
		T_STAGE2 m02 = T_STAGE1(a00*a21) - T_STAGE1(a20*a01);
		T_STAGE2 m12 = T_STAGE1(a10*a21) - T_STAGE1(a20*a11);
		
// 		std::cout << "m01=" << m01 << std::endl;
// 		std::cout << "m02=" << m02 << std::endl;
// 		std::cout << "m12=" << m12 << std::endl;
		
		//mij have 2n+1 bits
		
		// Now compute the minors of rank 3
		T_STAGE3 m012 = T_STAGE3(m01*a22) - T_STAGE3(m02*a12) + T_STAGE3(m12*a02);
		
// 		std::cout << "m012=" << m012 << std::endl;
		
		//mij*akl has 2n+1+n = 3n+1 bits
		//Thus m012 needs 3n+2 bits
		if (m012 < 0) {
			return CGAL::NEGATIVE;
		}
		else if (m012 > 0) {
			return CGAL::POSITIVE;
		}
		else {
			return CGAL::ZERO;
		}
	}
private:
	Is_in_auxiliary_triangle m_iat;
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

template<typename T_LINEAR_KERNEL>
struct EpsBasedAuxPointsTraits {
	using LinearKernel = T_LINEAR_KERNEL;
	using FT = typename LinearKernel::FT;
	using Point_3 = Point_sp<T_LINEAR_KERNEL>;
};

template<typename T_LINEAR_KERNEL>
struct AuxPointsGenerator {
public:
	using LinearKernel = T_LINEAR_KERNEL;
	using Point_3 = Point_sp<LinearKernel>;
	using FT = typename LinearKernel::FT;
// 	static constexpr uint32_t max_exponent = (static_cast<uint32_t>(1) << (Point_3::BitSizes::EXPONENT))-1;
	static constexpr uint32_t max_exponent = 60;
public:
	class Generate_auxiliary_point final {
	public:
		Generate_auxiliary_point() {}
	public:
		inline Point_3 operator()(AuxPointSelector s) const {
			Point_3 p;
			if (s == AuxPointSelector::LOWER) {
				p = Point_3(0,0,1, ratss::SP_DIM3_NEGATIVE);
			}
			else {
				p.set_pos(ratss::SP_DIM3_POSITIVE);
				p.set_exponent(max_exponent);
				switch(s) {
					case AuxPointSelector::UPPER_0:
						p.set_numerator0(0);
						p.set_numerator1(1);
						break;
					case AuxPointSelector::UPPER_1:
						p.set_numerator0(-1);
						p.set_numerator1(-1);
						break;
					case AuxPointSelector::UPPER_2:
						p.set_numerator0(1);
						p.set_numerator1(-1);
						break;
					default:
						throw std::runtime_error("AuxPointsGenerator::Generate_auxiliary_point: Invalid selector");
						break;
				}
			}
			assert(p == Point_3(p.point3()));
			return p;
		};
	};
	
	class Is_in_auxiliary_triangle final {
	public:
		Is_in_auxiliary_triangle() {}
	public:
		inline bool operator()(Point_3 const & p) const {
			return (p.pos() == ratss::SP_DIM3_POSITIVE &&
					p.exponent() >= max_exponent &&
					std::abs(p.numerator0()) <= 1 &&
					std::abs(p.numerator1()) <= 1
			);
		};
	};

	class Is_auxiliary_point final {
	public:
		Is_auxiliary_point() {}
	public:
		bool operator()(Point_3 const & p) const {
			return (p.pos() == ratss::SP_DIM3_NEGATIVE && p.numerator0() == 0 && p.numerator1() == 0) || m_iat(p);
		};
	private:
		Is_in_auxiliary_triangle m_iat;
	};
	
	//Is true if the given point is on the sphere and outside of the auxiliary triangle
	class Is_valid_point_on_sphere final {
	public:
		Is_valid_point_on_sphere(Is_auxiliary_point const & iap) : m_iap(iap) {}
	public:
		bool operator()(Point_3 const & p) const {
			return p.pos() != 0 && !m_iap(p);
		}
	private:
		Is_auxiliary_point m_iap;
	};
public:
	AuxPointsGenerator() {}
	virtual ~AuxPointsGenerator() {}
public:
	Is_auxiliary_point is_auxiliary_point_object() const {
		return Is_auxiliary_point();
	}
	
	Is_in_auxiliary_triangle is_in_auxiliary_triangle_object() const {
		return Is_in_auxiliary_triangle();
	}
	
	Is_valid_point_on_sphere is_valid_point_on_sphere_object() const {
		return Is_valid_point_on_sphere(is_auxiliary_point_object());
	}
	
	Generate_auxiliary_point generate_auxiliary_point_object() const {
		return Generate_auxiliary_point();
	}
};

} //end namespace detail::Kernel_sp 


template<
	typename T_LINEAR_KERNEL,
// 	typename T_AUX_POINT_GENERATOR=detail::EpsBasedAuxPoints<detail::Kernel_sp::EpsBasedAuxPointsTraits<T_LINEAR_KERNEL>>
	typename T_AUX_POINT_GENERATOR=detail::Kernel_sp::AuxPointsGenerator<T_LINEAR_KERNEL>
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
		return Orientation_2( MyBaseTrait::orientation_2_object(), MyBaseTrait::is_in_auxiliary_triangle_object() );
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
