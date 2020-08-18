#pragma once
#ifndef LIB_DTS2_CONSTRAINED_DELAUNAY_TRIANGULATION_WITH_EXACT_INTERSECTIONS_SPHERICAL_TRAITS_S2_H
#define LIB_DTS2_CONSTRAINED_DELAUNAY_TRIANGULATION_WITH_EXACT_INTERSECTIONS_SPHERICAL_TRAITS_S2_H

#include <libdts2/constants.h>

#include <CGAL/Exact_spherical_kernel_3.h>

#include <libratss/Conversion.h>
#include <libdts2/Constrained_delaunay_triangulation_base_traits_s2.h>

namespace LIB_DTS2_NAMESPACE {
namespace detail::spherical_traits {
	
class Is_auxiliary_point {
public:
	using MyKernel = CGAL::Exact_spherical_kernel_3;
	using MyLinearCDTKernel = Constrained_delaunay_triangulation_base_traits_s2<MyKernel::Linear_kernel>;
	using Point_2 = MyKernel::Circular_arc_point_3;
public:
	Is_auxiliary_point(MyLinearCDTKernel::Is_auxiliary_point const & base) : m_d(base) {}
public:
	bool operator()(Point_2 const & p) const {
		return m_d(p);
	}
private:
	MyLinearCDTKernel::Is_auxiliary_point m_d;
};

class Is_in_auxiliary_triangle {
public:
	using MyKernel = CGAL::Exact_spherical_kernel_3;
	using MyLinearCDTKernel = Constrained_delaunay_triangulation_base_traits_s2<MyKernel::Linear_kernel>;
	using Point_2 = MyKernel::Circular_arc_point_3;
public:
	Is_in_auxiliary_triangle(MyLinearCDTKernel::Is_in_auxiliary_triangle const & base) : m_d(base) {}
public:
	bool operator()(Point_2 const & p) const {
		return m_d(p);
	}
private:
	MyLinearCDTKernel:: Is_in_auxiliary_triangle m_d;
};

class Is_valid_point_on_sphere {
public:
	using MyKernel = CGAL::Exact_spherical_kernel_3;
	using MyLinearCDTKernel = Constrained_delaunay_triangulation_base_traits_s2<MyKernel::Linear_kernel>;
	using Point_2 = MyKernel::Circular_arc_point_3;
public:
	Is_valid_point_on_sphere(MyLinearCDTKernel::Is_auxiliary_point const & base) : m_iap(base) {}
public:
	bool operator()(Point_2 const & p) const {
		return !m_iap(p);
	}
private:
	MyLinearCDTKernel::Is_auxiliary_point m_iap;
};

class Generate_auxiliary_point {
public:
	using MyKernel = CGAL::Exact_spherical_kernel_3;
	using MyLinearCDTKernel = Constrained_delaunay_triangulation_base_traits_s2<MyKernel::Linear_kernel>;
	using Point_2 = MyKernel::Circular_arc_point_3;
public:
	Generate_auxiliary_point(MyLinearCDTKernel::Generate_auxiliary_point const & base) : m_d(base) {}
public:
	Point_2 operator()(AuxPointSelector s) const {
		return Point_2( m_d(s) );
	}
private:
	MyLinearCDTKernel:: Generate_auxiliary_point m_d;
};
	
class Triangle {
public:
	using K = CGAL::Exact_spherical_kernel_3;
	using Point_3 = K::Circular_arc_point_3;
public:
	Triangle(Point_3 const & a, Point_3 const & b, Point_3 const & c) : m_d({a, b, c}) {}
public:
	const Point_3 & vertex(int i) const { return m_d.at(i); }
private:
	std::array<Point_3, 3> m_d;
};

class Construct_triangle {
public:
	using K = CGAL::Exact_spherical_kernel_3;
	using Point_3 = K::Circular_arc_point_3;
public:
	Triangle operator()(Point_3 const & a, Point_3 const & b, Point_3 const & c) const {
		return Triangle(a, b, c);
	}
};

class Less_x: private CGAL::Exact_spherical_kernel_3::Compare_x_3 {
public:
	using MyBaseClass = CGAL::Exact_spherical_kernel_3::Compare_x_3;
public:
	template<typename T>
	bool operator()(T const & a, T const & b) const {
		return MyBaseClass::operator()(a, b) == CGAL::SMALLER;
	}
};

class Less_y: private CGAL::Exact_spherical_kernel_3::Compare_y_3 {
public:
	using MyBaseClass = CGAL::Exact_spherical_kernel_3::Compare_y_3;
public:
	template<typename T>
	bool operator()(T const & a, T const & b) const {
		return MyBaseClass::operator()(a, b) == CGAL::SMALLER;
	}
};

class Less_z: private CGAL::Exact_spherical_kernel_3::Compare_z_3 {
public:
	using MyBaseClass = CGAL::Exact_spherical_kernel_3::Compare_z_3;
public:
	template<typename T>
	bool operator()(T const & a, T const & b) const {
		return MyBaseClass::operator()(a, b) == CGAL::SMALLER;
	}
};

class Side_of_oriented_circle_s2 {
public:
	using K = CGAL::Exact_spherical_kernel_3;
	using LK = K::Linear_kernel;
	using Point_3 = LK::Point_3;
	using Point_s = K::Circular_arc_point_3;
public:
	CGAL::Orientation operator()(const Point_s & p, const Point_s & q, const Point_s & r, const Point_s & s) const {
		return ot3(p, q, r, s);
	}
private:
	CGAL::Orientation ot3(const Point_s & p, const Point_s & q, const Point_s & r, const Point_s & s) const {
		return CGAL::orientationC3<K::Root_of_2>(p.x(), p.y(), p.z(),
											  q.x(), q.y(), q.z(),
											  r.x(), r.y(), r.z(),
											  s.x(), s.y(), s.z()
								   );
	}
};

class Orientation_s2 {
public:
	using K = CGAL::Exact_spherical_kernel_3;
	using LK = K::Linear_kernel;
	using Point_3 = LK::Point_3;
	using Point_s = K::Circular_arc_point_3;
public:
	Orientation_s2(Is_in_auxiliary_triangle const & v) : m_iat(v) {}
public:
	CGAL::Orientation operator()(const Point_s & p, const Point_s & q, const Point_s & t) const {
		Point_s origin(LIB_DTS2_ORIGIN);
		CGAL::Sign oriented_side = ot3(p, q, origin, t);
		
		if (oriented_side != CGAL::Sign::COLLINEAR && m_iat(p) && m_iat(q)) {
			//the point always has to be on the opposite side of the infinite vertex
			oriented_side = - ot3(p, q, origin, Point_s(Point_3(0, 0, 1)));
		}
		return oriented_side;
	}
private:
	CGAL::Orientation ot3(const Point_s & p, const Point_s & q, const Point_s & r, const Point_s & s) const {
		return CGAL::orientationC3<K::Root_of_2>(p.x(), p.y(), p.z(),
											  q.x(), q.y(), q.z(),
											  r.x(), r.y(), r.z(),
											  s.x(), s.y(), s.z()
								   );
	}
private:
	Is_in_auxiliary_triangle m_iat;
};

class Construct_point {
public:
	using K = CGAL::Exact_spherical_kernel_3;
	using Point_3 = K::Circular_arc_point_3;
public:
	Point_3 const & operator()(Point_3 const & p) const { return p; }
	
	template<typename T>
	Point_3 operator()(T && p) const { return Point_3(std::forward<T>(p)); }
};

class Construct_segment_s2 {
public:
	using K = CGAL::Exact_spherical_kernel_3;
	using Point_3 = K::Circular_arc_point_3;
	using MyBase = K::Construct_circular_arc_3;
	using Segment = K::Circular_arc_3;
	using Circel_3 = K::Circle_3;
public:
	Construct_segment_s2() {}
public:
	Segment operator()(Point_3 const & src, Point_3 const & tgt) const {
		throw std::runtime_error("Unimplemented function");
		return Segment();
	}
private:
	MyBase m_d;
};

class Intersect_s2 {
public:
	using K = CGAL::Exact_spherical_kernel_3;
	using Point_3 = K::Circular_arc_point_3;
	using MyBase = K::Construct_circular_arc_3;
	using Segment = K::Circular_arc_3;
public:
	CGAL::Object operator()(Segment const & a, Segment const & b) const {
		throw std::runtime_error("Unimplemented function");
		return CGAL::Object();
	}
	
};
	
}//end namespace detail

class Constrained_delaunay_triangulation_with_exact_intersections_spherical_traits_s2 {
public:
	using MyKernel = CGAL::Exact_spherical_kernel_3;
	using MyLinearCDTKernel = Constrained_delaunay_triangulation_base_traits_s2<MyKernel::Linear_kernel>;
public:
	using FT = MyKernel::FT;
	using Oriented_side = MyKernel::Oriented_side;
	using Comparison_result = MyKernel::Comparison_result;
	using Orientation = MyKernel::Orientation;
	
	
	using Point_2 = MyKernel::Circular_arc_point_3;
	using Segment_2 = MyKernel::Circular_arc_3;
	using Triangle_2 = detail::spherical_traits::Triangle;
	using Construct_point_2 = detail::spherical_traits::Construct_point;
	using Construct_segment_2 = detail::spherical_traits::Construct_segment_s2;
	using Construct_triangle_2 = detail::spherical_traits::Construct_triangle;
	using Less_x_2 = detail::spherical_traits::Less_x;
	using Less_y_2 = detail::spherical_traits::Less_y;
	using Compare_x_2 = MyKernel::Compare_x_3;
	using Compare_y_2 = MyKernel::Compare_y_3;
	using Orientation_2 = detail::spherical_traits::Orientation_s2;
	using Side_of_oriented_circle_2 = detail::spherical_traits::Side_of_oriented_circle_s2;
	
	//the following is needed for spatial sorting
	using Less_x_3 = detail::spherical_traits::Less_x;
	using Less_y_3 = detail::spherical_traits::Less_y;
	using Less_z_3 = detail::spherical_traits::Less_z;
	
	using Point_3 = Point_2;
	using Point = Point_3;
	
	using Do_intersect_2 = MyKernel::Do_intersect_3;
	using Intersect_2 = detail::spherical_traits::Intersect_s2;

	//special stuff
	using Is_auxiliary_point = detail::spherical_traits::Is_auxiliary_point;
	using Is_in_auxiliary_triangle = detail::spherical_traits::Is_in_auxiliary_triangle;
	using Is_valid_point_on_sphere = detail::spherical_traits::Is_valid_point_on_sphere;
	using Generate_auxiliary_point = detail::spherical_traits::Generate_auxiliary_point;
	
public:
	class Project_on_sphere: private MyLinearCDTKernel::Project_on_sphere {
	public:
		using MyBaseClass = MyLinearCDTKernel::Project_on_sphere;
	private:
		template<typename T>
		using Conversion = LIB_RATSS_NAMESPACE::Conversion<T>;
	public:
		Project_on_sphere(const MyBaseClass & v) : MyBaseClass(v) {}
	public:
		using MyBaseClass::projector;
	public:
		Point_3 operator()(const Point_3 & /*v*/) const {
			throw std::runtime_error("Projecting a circular arc point is not supported");
			return Point_3();
		}
		Point_3 operator()(Point_3 & /*v*/) const {
			throw std::runtime_error("Projecting a circular arc point is not supported");
			return Point_3();
		}
		Point_3 operator()(FT const & x, FT const & y, FT const & z) const {
			auto sqLen = x*x + y*y + z*z;
			if (sqLen == 1) {
				return Point_3(MyLinearCDTKernel::Point_3(x, y, z));
			}
			
			mpq_class sqLenQ( Conversion<FT>::toMpq(sqLen) );
			std::size_t sqLenPrec = projector().calc().maxBitCount(sqLenQ);
			mpfr::mpreal sqLenF(Conversion<mpq_class>::toMpreal(sqLenQ, sqLenPrec));
			mpfr::mpreal lenF = projector().calc().sqrt(sqLenF);
			mpfr::mpreal xf(Conversion<FT>::toMpreal(x, sqLenPrec));
			mpfr::mpreal yf(Conversion<FT>::toMpreal(y, sqLenPrec));
			mpfr::mpreal zf(Conversion<FT>::toMpreal(z, sqLenPrec));
			mpq_class xq, yq, zq;
			projector().snap(xf, yf, zf, xq, yq, zq, significands());
			return Point_3( MyLinearCDTKernel::Point_3(
							Conversion<FT>::moveFrom(xq),
							Conversion<FT>::moveFrom(yq),
							Conversion<FT>::moveFrom(zq)
			));
		}
		template<typename... T>
		Point_3 operator()(T const & ... t) const {
			return Point_3( MyBaseClass::operator()(t...) );
		}
	};
public:
	Constrained_delaunay_triangulation_with_exact_intersections_spherical_traits_s2();
	Constrained_delaunay_triangulation_with_exact_intersections_spherical_traits_s2(MyLinearCDTKernel const & lk) : m_ltraits(lk) {}
	Constrained_delaunay_triangulation_with_exact_intersections_spherical_traits_s2(FT const & epsZ, int significands) : m_ltraits(epsZ, significands) {}
public:
	
	Side_of_oriented_circle_2 side_of_oriented_circle_2_object () const {
		return Side_of_oriented_circle_2();
	}
	
// 	Compare_distance_2 compare_distance_2_object () const {
// 		return baseTraits().compare_distance_3_object();
// 	}
	
	Orientation_2 orientation_2_object () const {
		return Orientation_2(is_in_auxiliary_triangle_object());
	}
	
	Compare_x_2 compare_x_2_object() const {
		return Compare_x_2();
	}
	
	Compare_y_2 compare_y_2_object() const {
		return Compare_y_2();
	}
	
	Less_x_3 less_x_3_object() const {
		return Less_x_3();
	}

	Less_y_3 less_y_3_object() const {
		return Less_y_3();
	}

	Less_z_3 less_z_3_object() const {
		return Less_z_3();
	}
	
	//stuff to support intersections
	Construct_segment_2 construct_segment_2_object() const {
		return Construct_segment_2();
	}

	Construct_point_2 construct_point_2_object() const {
		return Construct_point_2();
	}
	
	Less_x_2 less_x_2_object() const {
		return Less_x_2();
	}

	Less_y_2 less_y_2_object() const {
		return Less_y_2();
	}
	
	Do_intersect_2 do_intersect_2_object() const {
		return m_traits.do_intersect_3_object();
	}
	
	Intersect_2 intersect_2_object() const {
		return Intersect_2();
	}
	
	Project_on_sphere project_on_sphere_object() const {
		return Project_on_sphere(m_ltraits.project_on_sphere_object());
	}

	Project_on_sphere project_on_sphere_object(int _significands) const {
		return Project_on_sphere( m_ltraits.project_on_sphere_object(_significands) );
	}
	
	Is_auxiliary_point is_auxiliary_point_object() const {
		return  Is_auxiliary_point( m_ltraits.is_auxiliary_point_object() );
	}
	
	Is_in_auxiliary_triangle is_in_auxiliary_triangle_object() const {
		return Is_in_auxiliary_triangle( m_ltraits.is_in_auxiliary_triangle_object() );
	}
	
	Is_valid_point_on_sphere is_valid_point_on_sphere_object() const {
		return Is_valid_point_on_sphere( m_ltraits.is_auxiliary_point_object() );
	}
	
	Generate_auxiliary_point generate_auxiliary_point_object() const {
		return Generate_auxiliary_point( m_ltraits.generate_auxiliary_point_object() );
	}
private:
	MyKernel m_traits;
	MyLinearCDTKernel m_ltraits;
};


}//end namespace


#endif
