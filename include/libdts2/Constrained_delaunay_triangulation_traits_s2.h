#pragma once
#ifndef LIB_DTS2_CONSTRAINED_DELAUNAY_TRIANGULATION_TRAITS_S2_H
#define LIB_DTS2_CONSTRAINED_DELAUNAY_TRIANGULATION_TRAITS_S2_H

#include <libdts2/constants.h>
#include <libdts2/Constrained_delaunay_triangulation_base_traits_s2.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

namespace LIB_DTS2_NAMESPACE {

template<
	typename T_LINEAR_KERNEL,
	typename T_AUX_POINT_GENERATOR,
	typename T_POINT
>
class Constrained_delaunay_triangulation_traits_s2;

} //end namespace LIB_DTS2_NAMESPACE

namespace CGAL {
namespace internal {

	//this one is needed since Constrained_delaunay_triangulation_traits_s2 sets Point_2 = Point_3
	//and then we get a conflict of this function in the original code

    template <class ConcurrencyTag = Sequential_tag, class RandomAccessIterator, class Policy, typename T_LINEAR_KERNEL, typename T_AUX_POINT_GENERATOR, typename T_POINT>
    void spatial_sort (RandomAccessIterator begin, RandomAccessIterator end,
						const LIB_DTS2_NAMESPACE::Constrained_delaunay_triangulation_traits_s2<T_LINEAR_KERNEL, T_AUX_POINT_GENERATOR, T_POINT> & k, 
						Policy /*policy*/,
						std::ptrdiff_t threshold_hilbert,
						std::ptrdiff_t threshold_multiscale,
						double ratio);

} //end namespace internal
} //end namespace CGAL

namespace LIB_DTS2_NAMESPACE {

template<
	typename T_LINEAR_KERNEL,
	typename T_AUX_POINT_GENERATOR=detail::EpsBasedAuxPoints<T_LINEAR_KERNEL>,
	typename T_POINT = typename T_LINEAR_KERNEL::Point_3
>
class Constrained_delaunay_triangulation_traits_s2:
	public Constrained_delaunay_triangulation_base_traits_s2<T_LINEAR_KERNEL, T_AUX_POINT_GENERATOR, T_POINT>
{
public:
	using MyBaseTrait = Constrained_delaunay_triangulation_base_traits_s2<T_LINEAR_KERNEL, T_AUX_POINT_GENERATOR, T_POINT>;
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
	
	using Compare_distance_2 = typename MyBaseTrait::Compare_distance_3;
	
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
public:
	class Project_on_sphere: public MyBaseTrait::Project_on_sphere {
	public:
		using MyBaseClass = typename MyBaseTrait::Project_on_sphere;
	private:
		template<typename T>
		using Conversion = LIB_RATSS_NAMESPACE::Conversion<T>;
	public:
		Project_on_sphere(const MyBaseClass & v) : MyBaseClass(v) {}
	public:
		using MyBaseClass::operator();
		Point_3 operator()(const Point_3 & v) const {
			FT sqLen(v.x()*v.x() + v.y()*v.y() + v.z()*v.z());
			if (sqLen == 1) {
				return v;
			}
			
			mpq_class sqLenQ = Conversion<FT>::toMpq(sqLen);
			std::size_t sqLenPrec = MyBaseClass::projector().calc().maxBitCount(sqLenQ);
			mpfr::mpreal sqLenF(Conversion<mpq_class>::toMpreal(sqLenQ, sqLenPrec));
			mpfr::mpreal lenF = MyBaseClass::projector().calc().sqrt(sqLenF);
			mpfr::mpreal xf(Conversion<FT>::toMpreal(v.x(), sqLenPrec));
			mpfr::mpreal yf(Conversion<FT>::toMpreal(v.y(), sqLenPrec));
			mpfr::mpreal zf(Conversion<FT>::toMpreal(v.z(), sqLenPrec));
			mpq_class xq, yq, zq;
			MyBaseClass::projector().snap(xf, yf, zf, xq, yq, zq, MyBaseClass::significands());
			return Point_3( Conversion<FT>::moveFrom(xq),
							Conversion<FT>::moveFrom(yq),
							Conversion<FT>::moveFrom(zq)
			);
		}
	};
public:
	///This does not correctly initialize this trait!
	Constrained_delaunay_triangulation_traits_s2() {}
	///@param epsilon set the value of the z-coordinate above which no points should exist
	Constrained_delaunay_triangulation_traits_s2(const FT & _epsilon, int _significands) :
	MyBaseTrait(_epsilon, _significands)
	{}
	
	Constrained_delaunay_triangulation_traits_s2(const Constrained_delaunay_triangulation_traits_s2 & other) :
	MyBaseTrait(other)
	{}
	
	~Constrained_delaunay_triangulation_traits_s2() {}
public:
	inline Project_on_sphere project_on_sphere_object() const {
		return Project_on_sphere( MyBaseTrait::project_on_sphere_object() );
	}
};


}//end namespace LIB_DTS2_NAMESPACE


#include <CGAL/spatial_sort.h>

namespace CGAL {
namespace internal {

    template <class ConcurrencyTag, class RandomAccessIterator, class Policy, typename T_LINEAR_KERNEL, typename T_AUX_POINT_GENERATOR, typename T_POINT>
    void spatial_sort (RandomAccessIterator begin, RandomAccessIterator end,
						const LIB_DTS2_NAMESPACE::Constrained_delaunay_triangulation_traits_s2<T_LINEAR_KERNEL, T_AUX_POINT_GENERATOR, T_POINT>& k,
						Policy, std::ptrdiff_t threshold_hilbert, std::ptrdiff_t threshold_multiscale, double ratio)
    {
		using MyTraits = LIB_DTS2_NAMESPACE::Constrained_delaunay_triangulation_traits_s2<T_LINEAR_KERNEL, T_AUX_POINT_GENERATOR, T_POINT>;
		using Sort = Hilbert_sort_3<MyTraits, Policy>;
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
