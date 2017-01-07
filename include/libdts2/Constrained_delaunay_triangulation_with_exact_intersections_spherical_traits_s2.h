#ifndef LIB_DTS2_CONSTRAINED_DELAUNAY_TRIANGULATION_WITH_EXACT_INTERSECTIONS_SPHERICAL_TRAITS_S2_H
#define LIB_DTS2_CONSTRAINED_DELAUNAY_TRIANGULATION_WITH_EXACT_INTERSECTIONS_SPHERICAL_TRAITS_S2_H
#pragma once

#include <libdts2/constants.h>

#include <CGAL/Exact_spherical_kernel_3.h>

#include <libratss/Conversion.h>
#include <libdts2/Constrained_delaunay_triangulation_base_traits_s2.h>

namespace LIB_DTS2_NAMESPACE {

class Constrained_delaunay_triangulation_with_exact_intersections_spherical_traits_s2 {
public:
	using MyKernel = CGAL::Exact_spherical_kernel_3;
	using MyLinearCDTKernel = Constrained_delaunay_triangulation_base_traits_s2<MyKernel::Linear_kernel>;
public:
	using FT = MyKernel::FT;
	using Oriented_side = MyKernel::Oriented_side;
	using Comparison_result = MyKernel::Comparison_result;
	using Orientation = MyKernel::Orientation;
	
	using Segment_2 = MyKernel::Circular_arc_3;
	using Segment = Segment_2;
	
	using Point_3 = MyLinearCDTKernel::Point_3;
	using Point_2 = MyLinearCDTKernel::Point_2;
	using Point = MyLinearCDTKernel::Point;
	
	using Compare_distance_2 = MyLinearCDTKernel::Compare_distance_2;
	
	//the following is needed for spatial sorting
	using Less_x_3 = MyLinearCDTKernel::Less_x_3;
	using Less_y_3 = MyLinearCDTKernel::Less_y_3;
	using Less_z_3 = MyLinearCDTKernel::Less_z_3;
	
	//this is needed for spatial sorting during insertion of constraints
	//cgal uses Spatial_sort_traits_adapter_2 to sort inserted vertices
	//this in turn needs the following predicates
	//maybe we should just specialize Spatial_sort_traits_adapter_2 for this Trait
	
	using Less_x_2 = MyLinearCDTKernel::Less_x_2;
	using Less_y_2 = MyLinearCDTKernel::Less_y_2;
	
	using Construct_segment_2 = MyKernel::Construct_circular_arc_3;
public:
	class Project_on_sphere: public MyLinearCDTKernel::Project_on_sphere {
	public:
		using MyBaseClass = MyLinearCDTKernel::Project_on_sphere;
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
			
			mpq_class sqLenQ( Conversion<FT>::toMpq(sqLen) );
			std::size_t sqLenPrec = projector().calc().maxBitCount(sqLenQ);
			mpfr::mpreal sqLenF(Conversion<mpq_class>::toMpreal(sqLenQ, sqLenPrec));
			mpfr::mpreal lenF = projector().calc().sqrt(sqLenF);
			mpfr::mpreal xf(Conversion<FT>::toMpreal(v.x(), sqLenPrec));
			mpfr::mpreal yf(Conversion<FT>::toMpreal(v.y(), sqLenPrec));
			mpfr::mpreal zf(Conversion<FT>::toMpreal(v.z(), sqLenPrec));
			mpq_class xq, yq, zq;
			projector().snap(xf, yf, zf, xq, yq, zq, precision());
			return Point_3( Conversion<FT>::moveFrom(xq),
							Conversion<FT>::moveFrom(yq),
							Conversion<FT>::moveFrom(zq)
			);
		}
	};
public:
	Constrained_delaunay_triangulation_with_exact_intersections_spherical_traits_s2();
private:
	   MyKernel m_traits;
};


}//end namespace


#endif