#ifndef LIB_DTS2_CONSTRAINED_DELAUNAY_TRIANGULATION_WITH_INTERSECTIONS_BASE_TRAITS_S2_H
#define LIB_DTS2_CONSTRAINED_DELAUNAY_TRIANGULATION_WITH_INTERSECTIONS_BASE_TRAITS_S2_H
#pragma once

#include <libdts2/constants.h>
#include <libdts2/Constrained_delaunay_triangulation_base_traits_s2.h>

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
	using Plane_3 = typename MyBaseTrait::Plane_3;
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
	//pointOnSphere(const Vector_3 & v) -> Point_4 returning the point on the sphere where v would point to
	//if it had length 1
	template<typename T_BASE>
	class Intersect_2: private T_BASE {
	private:
		using MyBaseClass = T_BASE;
		using MyParent = Constrained_delaunay_triangulation_with_intersections_base_traits_s2<T_KERNEL>;
	public:
		Intersect_2(const Intersect_3 & it3, const Orientation_3 & ot3) : m_it3(it3), m_ot3 (ot3) {}
		Intersect_2(const Intersect_3 & it3, const Orientation_3 & ot3, const MyBaseClass & base) : MyBaseClass(base), m_it3(it3), m_ot3 (ot3) {}
		CGAL::Object operator()(const Segment & a, const Segment & b) const {
			Plane_3 aPlane(a, Point_3(0, 0, 0) );
			Plane_3 bPlane(b, Point_3(0, 0, 0) );
			
			auto rayObj = m_it3(aPlane, bPlane);
			Ray_3 ray = CGAL::object_cast<Ray_3>(rayObj);
			
			MyParent::print(ray.source());
			
			Vector_3 vec = ray.to_vector();
			
			Point_3 p1( MyBaseClass::pointOnSphere(vec) );
			Point_3 p2(pointReflect(p1));
			
			//since we're on the sphere there are two possible points on the sphere
			//TODO: We use the one which is inside the circumcircle of the smallest of the 4 possible triangle
			//We have to make sure that the sphere-point is on the opposite site
			Point_3 result;
			Point_3 tg1, tg2, tg3;
			getSmallestTriangle(a, b, tg1, tg2, tg3);
			Orientation ot_1 = m_ot3 (tg1, tg2, tg3, p1);
			Orientation ot_2 = m_ot3 (tg1, tg2, tg3, p2);
			Orientation ot_origin = m_ot3 (tg1, tg2, tg3, Point_3(0, 0, 0));
			if (ot_origin == CGAL::POSITIVE) {
				//origin is on the wrong side -> all ot-test flip their sign
				ot_1 = - ot_1;
				ot_2 = - ot_2;
			}
			//by definition on has to be on the one side an the other one on the other side
			//except for degenerate cases
			assert(ot_1 != ot_2);
			
			if (ot_1 == CGAL::POSITIVE) {
				result = std::move(p1);
			}
			else {
				result = std::move(p2);
			}
			
			CGAL::Object obj;
			obj.assign(result);
			return obj;
		}
	protected:
		Point_3 pointReflect(const Point_3 & p) const {
			return Point_3(p.x(), p.y(), p.z());
		}
		void getSmallestTriangle(const Segment & a, const Segment & b, Point_3 & tp1, Point_3 & tp2, Point_3 & tp3) const {
			tp1 = a.source();
			tp2 = a.target();
			tp3 = b.source();
		}
	private:
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