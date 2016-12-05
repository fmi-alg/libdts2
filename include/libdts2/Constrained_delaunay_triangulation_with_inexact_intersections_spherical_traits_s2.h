#ifndef LIB_DTS2_CONSTRAINED_DELAUNAY_TRIANGULATION_WITH_INEXACT_INTERSECTIONS_SPHERICAL_TRAITS_S2_H
#define LIB_DTS2_CONSTRAINED_DELAUNAY_TRIANGULATION_WITH_INEXACT_INTERSECTIONS_SPHERICAL_TRAITS_S2_H
#pragma once

#include <libdts2/constants.h>
#include <CGAL/Exact_spherical_kernel_3.h>


namespace LIB_DTS2_NAMESPACE {

class Constrained_delaunay_triangulation_with_inexact_intersections_spherical_traits_s2 {
public:
	using MyBaseTrait = CGAL::Exact_spherical_kernel_3;
protected:
	using Compare_distance_3 = MyBaseTrait::Linear_kernel::Compare_distance_3;
public:
	using FT = MyBaseTrait::FT;
	using Oriented_side = MyBaseTrait::Oriented_side;
	using Comparison_result = MyBaseTrait::Comparison_result;
	using Orientation = MyBaseTrait::Orientation;
	
	using Segment_2 = MyBaseTrait::Circular_arc_3;
	using Segment = Segment_2;
	
	//we want our triangulation to live in cartesian 3D-space
	//the CDT lives in 2D euclidean space so we have to redefine Point_2;
	using Point_3 = MyBaseTrait::Point_3;
	using Point_2 = Point_3;
	using Point = Point_2;
	
	using Compare_distance_2 = Compare_distance_3;
	
	//the following is needed for spatial sorting
	using Less_x_3 = MyBaseTrait::Linear_kernel::Less_x_3;
	using Less_y_3 = MyBaseTrait::Linear_kernel::Less_y_3;
	using Less_z_3 = MyBaseTrait::Linear_kernel::Less_z_3;
	
	//this is needed for spatial sorting during insertion of constraints
	//cgal uses Spatial_sort_traits_adapter_2 to sort inserted vertices
	//this in turn needs the following predicates
	//maybe we should just specialize Spatial_sort_traits_adapter_2 for this Trait
	
	using Less_x_2 = Less_x_3;
	using Less_y_2 = Less_y_3;
	
	using Construct_segment_2 = MyBaseTrait::Construct_circular_arc_3;
public:
	Constrained_delaunay_triangulation_with_inexact_intersections_spherical_traits_s2();
private:
	
private:
	MyBaseTrait m_traits;
};

}//end namespace


#endif