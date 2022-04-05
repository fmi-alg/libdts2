#pragma once
#ifndef LIB_DTS2_CONSTRAINED_DELAUNAY_TRIANGULATION_WITH_EXACT_INTERSECTIONS_TRAITS_S2_H
#define LIB_DTS2_CONSTRAINED_DELAUNAY_TRIANGULATION_WITH_EXACT_INTERSECTIONS_TRAITS_S2_H

#include <libdts2/constants.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>

#include <libdts2/Constrained_delaunay_triangulation_with_intersections_base_traits_s2.h>

namespace LIB_DTS2_NAMESPACE {

class Constrained_delaunay_triangulation_with_exact_intersections_traits_s2;

} //end namespace LIB_DTS2_NAMESPACE

namespace CGAL {
namespace internal {

	//this one is needed since Constrained_delaunay_triangulation_traits_s2 sets Point_2 = Point_3
	//and then we get a conflict of this function in the original code
	
    template <class RandomAccessIterator, class Policy>
    void spatial_sort (
                       RandomAccessIterator begin, RandomAccessIterator end,
                       const LIB_DTS2_NAMESPACE::Constrained_delaunay_triangulation_with_exact_intersections_traits_s2 & k, 
		       Policy /*policy*/,
		       std::ptrdiff_t threshold_hilbert,
		       std::ptrdiff_t threshold_multiscale,
		       double ratio);

} //end namespace internal

template <class RandomAccessIterator, class Policy>
void spatial_sort (RandomAccessIterator begin, RandomAccessIterator end,
                   const LIB_DTS2_NAMESPACE::Constrained_delaunay_triangulation_with_exact_intersections_traits_s2 &k,
		   Policy policy,
		   std::ptrdiff_t threshold_hilbert=0,
		   std::ptrdiff_t threshold_multiscale=0,
		   double ratio=0.0)
{
    internal::spatial_sort(begin, end, k, policy, threshold_hilbert,threshold_multiscale,ratio);
}

} //end namespace CGAL


namespace LIB_DTS2_NAMESPACE {


///This traits class supports exact intersections of constrained segments
class Constrained_delaunay_triangulation_with_exact_intersections_traits_s2:
	public Constrained_delaunay_triangulation_with_intersections_base_traits_s2<
		CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt
	>
{
public:
	using MyBaseTrait = Constrained_delaunay_triangulation_with_intersections_base_traits_s2<CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt>;
public:
	using FT = MyBaseTrait::FT;
	using Oriented_side = MyBaseTrait::Oriented_side;
	using Comparison_result = MyBaseTrait::Comparison_result;
	using Orientation = MyBaseTrait::Orientation;
	
	using Segment_2 = MyBaseTrait::Segment_2;
	using Segment = MyBaseTrait::Segment;
	
	using Point_3 = MyBaseTrait::Point_3;
	using Point_2 = MyBaseTrait::Point_2;
	using Point = MyBaseTrait::Point;
	
	using Compare_distance_2 = Compare_distance_3;
	
	//the following is needed for spatial sorting
	using Less_x_3 = MyBaseTrait::Less_x_3;
	using Less_y_3 = MyBaseTrait::Less_y_3;
	using Less_z_3 = MyBaseTrait::Less_z_3;
	
	//this is needed for spatial sorting during insertion of constraints
	//cgal uses Spatial_sort_traits_adapter_2 to sort inserted vertices
	//this in turn needs the following predicates
	//maybe we should just specialize Spatial_sort_traits_adapter_2 for this Trait
	
	using Less_x_2 = MyBaseTrait::Less_x_2;
	using Less_y_2 = MyBaseTrait::Less_y_2;
	
	using Construct_segment_2 = MyBaseTrait::Construct_segment_2;
public:
	class Project_on_sphere: public MyBaseTrait::Project_on_sphere {
	public:
		using MyBaseClass = MyBaseTrait::Project_on_sphere;
	public:
		Project_on_sphere(const MyBaseClass & base) : MyBaseClass(base) {}
	public:
		using MyBaseClass::operator();
		Point_3 operator()(const Point_3 & v) const {
			if ( CGAL::compare_squared_distance(v, Point_3(0,0,0), FT(1)) != CGAL::Comparison_result::EQUAL) {
				FT len( CGAL::sqrt(v.x()*v.x() + v.y()*v.y() + v.z()*v.z()) );
				return Point_3(v.x() / len, v.y() / len, v.z() / len);
			}
			else {
				return v;
			}
		}
	};
protected: //own implementations not support by the base traits
	class Intersect_base_2 {
	public:
		Intersect_base_2() {}
		Intersect_base_2(const Intersect_base_2 &) {}
		~Intersect_base_2() {}
	protected:
		Point_3 pointOnSphere(const Vector_3 & v) const {
			FT len( CGAL::sqrt(v.x()*v.x() + v.y()*v.y() + v.z()*v.z()) );
			return Point_3(v.x() / len, v.y() / len, v.z() / len);
		}
	};
public:
	using Intersect_2 = MyBaseTrait::Intersect_2<Intersect_base_2>;
public:
	///This does not correctly initialize this trait!
	Constrained_delaunay_triangulation_with_exact_intersections_traits_s2() {}
	///@param epsilon set the value of the z-coordinate above which no points should exist
	Constrained_delaunay_triangulation_with_exact_intersections_traits_s2(const FT & _epsilon, int _significands) :
	MyBaseTrait(_epsilon, _significands)
	{}
	
	Constrained_delaunay_triangulation_with_exact_intersections_traits_s2(const Constrained_delaunay_triangulation_with_exact_intersections_traits_s2 & other) :
	MyBaseTrait(other)
	{}

	~Constrained_delaunay_triangulation_with_exact_intersections_traits_s2() {}
public: //predicate creation functions
	inline Intersect_2 intersect_2_object() const {
		return Intersect_2(MyBaseTrait::do_intersect_3_object(), MyBaseTrait::intersect_3_object(), MyBaseTrait::orientation_3_object());
	}
	
	inline Project_on_sphere project_on_sphere_object() const {
		return Project_on_sphere(MyBaseTrait::project_on_sphere_object());
	}
};

}//end LIB_DTS2_NAMESPACE

#include <CGAL/spatial_sort.h>

namespace CGAL {
namespace internal {

    template <class RandomAccessIterator, class Policy>
    void spatial_sort (RandomAccessIterator begin, RandomAccessIterator end,
						const LIB_DTS2_NAMESPACE::Constrained_delaunay_triangulation_with_exact_intersections_traits_s2 & k,
						Policy, std::ptrdiff_t threshold_hilbert, std::ptrdiff_t threshold_multiscale, double ratio)
    {
      typedef Hilbert_sort_3<LIB_DTS2_NAMESPACE::Constrained_delaunay_triangulation_with_exact_intersections_traits_s2, Policy> Sort;
        boost::rand48 random;
        boost::random_number_generator<boost::rand48> rng(random);
#if defined(CGAL_HILBERT_SORT_WITH_MEDIAN_POLICY_CROSS_PLATFORM_BEHAVIOR)
        CGAL::random_shuffle(begin,end, rng);
#else
        std::random_shuffle(begin,end, rng);
#endif

	if (threshold_hilbert==0) threshold_hilbert=8;
	if (threshold_multiscale==0) threshold_multiscale=64;
	if (ratio==0.0) ratio=0.125;

        (Multiscale_sort<Sort> (Sort (k, threshold_hilbert), 
				threshold_multiscale, ratio)) (begin, end);
    }

}} //end namespace CGAL::internal

#endif
