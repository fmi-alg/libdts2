#pragma once
#ifndef LIB_DTS2_CONSTRAINED_DELAUNAY_TRIANGULATION_WITH_INEXACT_INTERSECTIONS_TRAITS_S2_H
#define LIB_DTS2_CONSTRAINED_DELAUNAY_TRIANGULATION_WITH_INEXACT_INTERSECTIONS_TRAITS_S2_H

#include <libdts2/constants.h>
#include <libdts2/Constrained_delaunay_triangulation_with_intersections_base_traits_s2.h>
#include <libdts2/Kernel_sp/Point_sp.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <libratss/CGAL/ExtendedInt64Cartesian.h>


namespace LIB_DTS2_NAMESPACE {

template<
	typename T_LINEAR_KERNEL,
	typename T_AUX_POINT_GENERATOR,
	typename T_POINT
>
class Constrained_delaunay_triangulation_with_inexact_intersections_base_traits_s2;

} //end namespace LIB_DTS2_NAMESPACE

namespace CGAL {
namespace internal {

	//this one is needed since Constrained_delaunay_triangulation_traits_s2 sets Point_2 = Point_3
	//and then we get a conflict of this function in the original code
	
    template <class RandomAccessIterator, class Policy, typename Kernel, typename T_AUX_POINT_GENERATOR, typename T_POINT>
    void spatial_sort (
                       RandomAccessIterator begin, RandomAccessIterator end,
                       const LIB_DTS2_NAMESPACE::Constrained_delaunay_triangulation_with_inexact_intersections_base_traits_s2<Kernel, T_AUX_POINT_GENERATOR, T_POINT> & k, 
		       Policy /*policy*/,
		       std::ptrdiff_t threshold_hilbert,
		       std::ptrdiff_t threshold_multiscale,
		       double ratio);

} //end namespace internal

} //end namespace CGAL


namespace LIB_DTS2_NAMESPACE {

///This traits class supports exact intersections of constrained segments
template<
	typename T_LINEAR_KERNEL,
	typename T_AUX_POINT_GENERATOR=detail::EpsBasedAuxPoints<T_LINEAR_KERNEL>,
	typename T_POINT = typename T_LINEAR_KERNEL::Point_3
>
class Constrained_delaunay_triangulation_with_inexact_intersections_base_traits_s2:
	public Constrained_delaunay_triangulation_with_intersections_base_traits_s2<T_LINEAR_KERNEL, T_AUX_POINT_GENERATOR, T_POINT>
{
public:
	using MyBaseTrait = Constrained_delaunay_triangulation_with_intersections_base_traits_s2<T_LINEAR_KERNEL, T_AUX_POINT_GENERATOR, T_POINT>;
	using MyKernel = typename MyBaseTrait::MyKernel;
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
	
	using Vector_3 = typename MyBaseTrait::Vector_3;
	
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
	
	//the following is needed for inexact intersection calculation
	using Compute_squared_distance_2 = typename MyKernel::Compute_squared_distance_3;
	using Construct_line_2 = typename MyKernel::Construct_line_3;
	using Line_2 = typename MyKernel::Line_3;
public:
	class Project_on_sphere: public MyBaseTrait::Project_on_sphere {
	public:
		using MyBaseClass = typename MyBaseTrait::Project_on_sphere;
	private:
		template<typename T>
		using Conversion = LIB_RATSS_NAMESPACE::Conversion<T>;
	public:
		Project_on_sphere(const MyBaseClass & base) : MyBaseClass(base) {}
	public:
		using MyBaseClass::operator();
		Point_3 operator()(const Point_3 & v) const {
			return (*this)(v.x(), v.y(), v.z());
		}
		
		//This variant is needed for Point_sp since we can only store points on the sphere
		Point_3 operator()(FT const & x, FT const & y, FT const & z) const {
			FT sqLen(x*x + y*y + z*z);
			if (sqLen == 1) {
				return Point_3(x, y, z);
			}
			mpq_class sqLenQ( Conversion<FT>::toMpq(sqLen) );
			
			{
				mpz_class tmp = sqrt(sqLenQ.get_num());
				if (tmp*tmp == sqLenQ.get_num()) {
					mpz_class tmp2 = sqrt(sqLenQ.get_den());
					if (tmp2*tmp2 == sqLenQ.get_den()) {
						FT len = Conversion<FT>::moveFrom(mpq_class(tmp)/mpq_class(tmp2));
						return Point_3( x/len,
										y/len,
										z/len
						);
					}
				}
			}
			
			int calcPrec = std::max<std::size_t>(MyBaseClass::projector().calc().maxBitCount(sqLenQ), MyBaseClass::significands()*2);
			
			mpfr::mpreal sqLenF(Conversion<mpq_class>::toMpreal(sqLenQ, calcPrec));
			mpfr::mpreal lenF = MyBaseClass::projector().calc().sqrt(sqLenF);
			
			mpfr::mpreal xf(Conversion<FT>::toMpreal(x, calcPrec));
			mpfr::mpreal yf(Conversion<FT>::toMpreal(y, calcPrec));
			mpfr::mpreal zf(Conversion<FT>::toMpreal(z, calcPrec));
			
			xf /= lenF;
			yf /= lenF;
			zf /= lenF;
			
			mpq_class xq, yq, zq;
			MyBaseClass::projector().snap(xf, yf, zf, xq, yq, zq, MyBaseClass::significands(), LIB_RATSS_NAMESPACE::ST_FX | LIB_RATSS_NAMESPACE::ST_PLANE);
			return Point_3( Conversion<FT>::moveFrom(xq),
							Conversion<FT>::moveFrom(yq),
							Conversion<FT>::moveFrom(zq)
			);
		}
	};
	
protected: //own implementations not support by the base traits
	class Intersect_base_2 {
	public:
		Intersect_base_2(const Project_on_sphere & pos) : m_pos(pos) {}
		Intersect_base_2(const Intersect_base_2 & other) : m_pos(other.m_pos) {}
		~Intersect_base_2() {}
	protected:
		inline Point_3 pointOnSphere(const Vector_3 & v) const {
			return m_pos(v.x(), v.y(), v.z());
		}
	private:
		Project_on_sphere m_pos;
	};
public:
	using Intersect_2 = typename MyBaseTrait::template Intersect_2<Intersect_base_2>;
protected:
	inline Intersect_base_2 intersect_base_2_object() const {
		return Intersect_base_2( project_on_sphere_object( intersectSignificands() ) );
	}
public:
	///This does not correctly initialize this trait!
	Constrained_delaunay_triangulation_with_inexact_intersections_base_traits_s2() : m_intersectSignificands(-1) {}
	///@param epsilon set the value of the z-coordinate above which no points should exist
	Constrained_delaunay_triangulation_with_inexact_intersections_base_traits_s2(const FT & _epsilon, int _significands, int _intersectSignificands = -1) :
	MyBaseTrait(_epsilon, _significands ),
	m_intersectSignificands(_intersectSignificands > 0 ? _intersectSignificands : _significands)
	{}
	
	Constrained_delaunay_triangulation_with_inexact_intersections_base_traits_s2(const Constrained_delaunay_triangulation_with_inexact_intersections_base_traits_s2 & other) :
	MyBaseTrait(other),
	m_intersectSignificands(other.m_intersectSignificands)
	{}

	~Constrained_delaunay_triangulation_with_inexact_intersections_base_traits_s2() {}
public: //predicate creation Constrained_delaunay_triangulation_with_inexact_intersections_traits_s2
	inline Intersect_2 intersect_2_object() const {
		return Intersect_2(MyBaseTrait::do_intersect_3_object(), MyBaseTrait::intersect_3_object(), MyBaseTrait::orientation_3_object(), intersect_base_2_object());
	}
	
	inline Construct_line_2 construct_line_2_object() const {
		return Construct_line_2();
	}
	
	inline Compute_squared_distance_2 compute_squared_distance_2_object() const {
		return Compute_squared_distance_2();
	}
	
	inline Project_on_sphere project_on_sphere_object() const {
		return Project_on_sphere(MyBaseTrait::project_on_sphere_object());
	}
	
	inline Project_on_sphere project_on_sphere_object(int _significands) const {
		return Project_on_sphere(MyBaseTrait::project_on_sphere_object(_significands));
	}
public:
	inline int intersectSignificands() const {
		return m_intersectSignificands;
	}
private:
	int m_intersectSignificands;
};

using Constrained_delaunay_triangulation_with_inexact_intersections_traits_s2 = 
	Constrained_delaunay_triangulation_with_inexact_intersections_base_traits_s2<CGAL::Exact_predicates_exact_constructions_kernel>;

///This trait uses ExtendedInt64q as number type which should in theory use less space if you snap with 31 Bits
using Constrained_delaunay_triangulation_with_inexact_intersections_traits_s2_64 =
// 	Constrained_delaunay_triangulation_with_inexact_intersections_base_traits_s2<CGAL::Simple_cartesian_extended_integer_kernel>;
	Constrained_delaunay_triangulation_with_inexact_intersections_base_traits_s2<CGAL::Filtered_simple_cartesian_extended_integer_kernel>;
// 	Constrained_delaunay_triangulation_with_inexact_intersections_base_traits_s2<CGAL::Filtered_lazy_cartesian_extended_integer_kernel>;
// 	Constrained_delaunay_triangulation_with_inexact_intersections_base_traits_s2<CGAL::Filtered_lazy_cartesian_extended_1024_integer_kernel>;
// 	Constrained_delaunay_triangulation_with_inexact_intersections_base_traits_s2<CGAL::Simple_cartesian_extended_1024_integer_kernel>;
// 	Constrained_delaunay_triangulation_with_inexact_intersections_base_traits_s2<CGAL::Exact_predicates_exact_constructions_extended_integer_kernel>;
	
using Constrained_delaunay_triangulation_with_inexact_intersections_traits_s2_sp =
	Constrained_delaunay_triangulation_with_inexact_intersections_base_traits_s2<
		CGAL::Filtered_simple_cartesian_extended_integer_kernel,
		detail::EpsBasedAuxPoints<CGAL::Filtered_simple_cartesian_extended_integer_kernel>,
		Point_sp<CGAL::Filtered_simple_cartesian_extended_integer_kernel>
	>;

}//end LIB_DTS2_NAMESPACE

#include <CGAL/spatial_sort.h>

namespace CGAL {
namespace internal {

    template <class RandomAccessIterator, class Policy, typename Kernel, typename T_AUX_POINT_GENERATOR, typename T_POINT>
    void spatial_sort (RandomAccessIterator begin, RandomAccessIterator end,
						const LIB_DTS2_NAMESPACE::Constrained_delaunay_triangulation_with_inexact_intersections_base_traits_s2<Kernel, T_AUX_POINT_GENERATOR, T_POINT> & k,
						Policy, std::ptrdiff_t threshold_hilbert, std::ptrdiff_t threshold_multiscale, double ratio)
    {
		using MyTraits = LIB_DTS2_NAMESPACE::Constrained_delaunay_triangulation_with_inexact_intersections_base_traits_s2<Kernel, T_AUX_POINT_GENERATOR, T_POINT>;
		
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
