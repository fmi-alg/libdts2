#ifndef LIB_DTS2_CONSTRAINED_DELAUNAY_TRIANGULATION_WITH_INTERSECTIONS_BASE_TRAITS_S2_H
#define LIB_DTS2_CONSTRAINED_DELAUNAY_TRIANGULATION_WITH_INTERSECTIONS_BASE_TRAITS_S2_H
#pragma once

#include <libdts2/constants.h>
#include <libdts2/Constrained_delaunay_triangulation_base_traits_s2.h>
#include <CGAL/intersections.h>

namespace LIB_DTS2_NAMESPACE {

template<typename T_KERNEL>
class Constrained_delaunay_triangulation_with_intersections_base_traits_s2:
	public Constrained_delaunay_triangulation_base_traits_s2<T_KERNEL>
{
public:
	using MyBaseTrait = Constrained_delaunay_triangulation_base_traits_s2<T_KERNEL>;
	using MyKernel = typename MyBaseTrait::MyKernel;
protected:
	using Segment_3 = typename MyBaseTrait::Segment_3;
	using Orientation_3 = typename MyBaseTrait::Orientation_3;

	using Compare_distance_3 = typename MyBaseTrait::Compare_distance_3;
	
	//needed by Compare_x_2 and Compare_y_2
	using Compare_x_3 = typename MyBaseTrait::Compare_x_3;
	using Compare_y_3 = typename MyBaseTrait::Compare_y_3;
	using Compare_z_3 = typename MyBaseTrait::Compare_z_3;
	
	
	using Triangle_3 = typename MyBaseTrait::Triangle_3;
	
	//stuff for intersection support
	using Construct_segment_3 = typename MyBaseTrait::Construct_segment_3;
	using Intersect_3 = typename MyBaseTrait::Intersect_3;
	using Do_intersect_3 = typename MyBaseTrait::Do_intersect_3;
	using Plane_3 = typename MyBaseTrait::Plane_3;
	using Line_3 = typename MyBaseTrait::Line_3;
	using Ray_3 = typename MyBaseTrait::Ray_3;
	using Vector_3 = typename MyBaseTrait::Vector_3;
public:
	using FT = typename MyBaseTrait::FT;
	using Oriented_side = typename MyBaseTrait::Oriented_side;
	using Comparison_result = typename MyBaseTrait::Comparison_result;
	using Orientation = typename MyBaseTrait::Orientation;
	
	using Segment_2 = typename MyBaseTrait::Segment_2;
	using Segment = typename MyBaseTrait::Segment;
	
	using Point_3 = typename MyBaseTrait::Point_3;
	using Point_2 = typename MyBaseTrait::Point_2;
	using Point = typename MyBaseTrait::Point;
	
	using Compare_distance_2 = Compare_distance_3;
	
	//the following is needed for spatial sorting
	using Less_x_3 = typename MyBaseTrait::Less_x_3;
	using Less_y_3 = typename MyBaseTrait::Less_y_3;
	using Less_z_3 = typename MyBaseTrait::Less_z_3;
	
	//this is needed for spatial sorting during insertion of constraints
	//cgal uses Spatial_sort_traits_adapter_2 to sort inserted vertices
	//this in turn needs the following predicates
	//maybe we should just specialize Spatial_sort_traits_adapter_2 for this Trait
	
	using Less_x_2 = typename MyBaseTrait::Less_x_2;
	using Less_y_2 = typename MyBaseTrait::Less_y_2;
	
	using Construct_segment_2 = typename MyBaseTrait::Construct_segment_2;
protected:
	using MyBaseTrait::print;
protected: //own implementations not support by the base traits
	//stuff needed for intersection, you have to provide a matching base class with a function
	//pointOnSphere(const Vector_3 & v) -> Point_3 returning the point on the sphere where v would point to
	//if it had length 1
	template<typename T_BASE>
	class Intersect_2: private T_BASE {
	private:
		using MyBaseClass = T_BASE;
		using MyParent = Constrained_delaunay_triangulation_with_intersections_base_traits_s2<T_KERNEL>;
		using Conversion = LIB_RATSS_NAMESPACE::Conversion<FT>;
	public:
		Intersect_2(const Do_intersect_3 & dit3, const Intersect_3 & it3, const Orientation_3 & ot3) :
		m_dit3(dit3), m_it3(it3), m_ot3 (ot3)
		{}
		Intersect_2(const Do_intersect_3 & dit3, const Intersect_3 & it3, const Orientation_3 & ot3, const MyBaseClass & base) :
		MyBaseClass(base), m_dit3(dit3), m_it3(it3), m_ot3 (ot3)
		{}
		CGAL::Object operator()(const Segment & a, const Segment & b) const {
// 			std::cerr << "Intersect_2 called -- BEGIN" << std::endl;
// 			std::cerr << "a=" << a << std::endl;
// 			std::cerr << "b=" << b << std::endl;
		
			Plane_3 aPlane(a, Point_3(0, 0, 0) );
			Plane_3 bPlane(b, Point_3(0, 0, 0) );
// 			std::cerr << "aPlane=" << aPlane << std::endl;
// 			std::cerr << "bPlane=" << aPlane << std::endl;
			
			auto xRes = m_it3(aPlane, bPlane);
			if (!boost::get<Line_3>(&*xRes)) {
				throw std::runtime_error("Trying to intersect segments on the same great circle is currenty unsupported");
			}
			
			const Line_3 * line3 = boost::get<Line_3>(&*xRes);
			
			//we need to check if the line passes through the segment a and b
			//This has to be done since the line3 represents the intersection points of the two great circles induces by segments a and b
			//Their intersection is not necessarily an intersection of the segments.
			if (!m_dit3(*line3, a) || !m_dit3(*line3, b)) {
				assert(false);
				return CGAL::Object();
			}
			
// 			std::cerr << "*line3=" << *line3 << std::endl;
			
			
			Vector_3 vec = line3->to_vector();
// 			std::cerr << "vec=" << vec << std::endl;
			
			
			Point_3 p1( MyBaseClass::pointOnSphere(vec) );
			Point_3 p2(pointReflect(p1));
			
// 			std::cerr << "p1=" << p1 << std::endl;
// 			std::cerr << "p2=" << p2 << std::endl;
			
			FT ds1 = minSqDist(a, b, p1);
			FT ds2 = minSqDist(a, b, p2);

// 			std::cerr << "ds1=" << ds1 << std::endl;
// 			std::cerr << "ds2=" << ds2 << std::endl;

			Point_3 result;
			if (ds1 < ds2) {
				result = std::move(p1);
			}
			else {
				result = std::move(p2);
			}
// 			assert(result != a.source() || line3->has_on(result));
// 			assert(result != b.source() || line3->has_on(result));
// 			assert(result != a.target() || line3->has_on(result));
// 			assert(result != b.target() || line3->has_on(result));
			
// 			std::cerr << "result=" << result << std::endl;
// 			std::cerr << "Intersect_2 called -- END" << std::endl;
// 			bool resultIsEndpoint = result == a.source() || result == b.source() || result == a.target() || result == b.target();
			std::cerr << "Intersect_2([";
			MyBaseTrait::print(std::cerr, a.source());
			std::cerr << ";";
			MyBaseTrait::print(std::cerr, a.target());
			std::cerr << "],[";
			MyBaseTrait::print(std::cerr, b.source());
			std::cerr << ";";
			MyBaseTrait::print(std::cerr, b.target());
			std::cerr  << "]) -> ";
			MyBaseTrait::print(std::cerr, result);
			std::cerr << std::endl;
			return CGAL::make_object(result);
		}
	protected:
		Point_3 pointReflect(const Point_3 & p) const {
			return Point_3(-p.x(), -p.y(), -p.z());
		}
		FT sq(const FT & ft) const {
			return ft*ft;
		}
		FT sqDist(const Point_3 & a, const Point_3 & b) const {
			return sq(a.x()-b.x()) + sq(a.y()-b.y()) + sq(a.z()-b.z());
		}
		FT minSqDist(const Segment & a, const Segment & b, const Point_3 & p) const {
			auto ds1 = sqDist(a.source(), p);
			auto ds2 = sqDist(a.target(), p);
			auto ds3 = sqDist(b.source(), p);
			auto ds4 = sqDist(b.target(), p);
			using std::min;
			return min(min(ds1, ds2), min(ds3, ds4));
		}
		
		//TODO:need implementation
		void getSmallestTriangle(const Segment & a, const Segment & b, Point_3 & tp1, Point_3 & tp2, Point_3 & tp3) const {
			tp1 = a.source();
			tp2 = a.target();
			tp3 = b.source();
		}
	private:
		Do_intersect_3 m_dit3;
		Intersect_3 m_it3;
		Orientation_3 m_ot3;
	};
public:
	///This does not correctly initialize this trait!
	Constrained_delaunay_triangulation_with_intersections_base_traits_s2() {}
	///@param epsilon set the value of the z-coordinate above which no points should exist
	Constrained_delaunay_triangulation_with_intersections_base_traits_s2(const FT & _epsilon, int _precision) :
	MyBaseTrait(_epsilon, _precision)
	{}
	
	Constrained_delaunay_triangulation_with_intersections_base_traits_s2(const Constrained_delaunay_triangulation_with_intersections_base_traits_s2 & other) :
	MyBaseTrait(other)
	{}

	~Constrained_delaunay_triangulation_with_intersections_base_traits_s2() {}
};

}//end LIB_DTS2_NAMESPACE

#endif