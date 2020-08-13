#pragma once
#ifndef LIB_DTS2_CONSTRAINED_DELAUNAY_TRIANGULATION_S2_H
#define LIB_DTS2_CONSTRAINED_DELAUNAY_TRIANGULATION_S2_H

#include <libdts2/constants.h>
#include <libdts2/Delaunay_triangulation_s2.h>
#include <libdts2/Constrained_Triangulation_base_s2.h>

//we have to include these first (due to some overloaded internal functions)
#include <libdts2/Constrained_delaunay_triangulation_with_exact_intersections_spherical_traits_s2.h>
#include <libdts2/Constrained_delaunay_triangulation_with_exact_intersections_traits_s2.h>
#include <libdts2/Constrained_delaunay_triangulation_with_inexact_intersections_traits_s2.h>
#include <libdts2/Kernel_sp/Kernel_sp.h>
#include <libdts2/Constrained_delaunay_triangulation_traits_s2.h>

#include <CGAL/Constrained_Delaunay_triangulation_2.h>

namespace LIB_DTS2_NAMESPACE {
namespace internal {
	template<typename T_KERNEL, typename T_VERTEX_INFO=void, typename T_FACE_INFO=void>
	class ConstrainedTriangulationDataStructureSelector {
	private:
		using Vb = typename VertexBaseSelector<T_KERNEL, T_VERTEX_INFO>::type;
		using Fb = typename FaceBaseSelector<T_KERNEL, T_FACE_INFO>::type;
		using CFb = CGAL::Constrained_triangulation_face_base_2<T_KERNEL, Fb>;
	public:
		using type = CGAL::Triangulation_data_structure_2<Vb, CFb>;
	};
} //end namespace internal


/** This is a constrained delaunay triangulation living on the sphere S2
	It uses CGAL::Constrained_Delaunay_triangulation_2 as bases
	Depending on the interstion tag the GeometricTraits need to obey several concepts
	CGAL::No_intersection_tag -> Constrained_delaunay_triangulation_traits_s2 concept
	CGAL::Exact_predicates_tag -> Constrained_delaunay_triangulation_with_inexact_intersections_traits_s2
	CGAL::Exact_intersections_tag -> Constrained_delaunay_triangulation_with_exact_intersections_traits_s2
**/
template<typename T_GEOM_TRAITS, typename T_TDS, typename T_ITAG>
class Constrained_delaunay_triangulation_s2:
	public Constrained_triangulation_base_s2<
		CGAL::Constrained_Delaunay_triangulation_2<T_GEOM_TRAITS, T_TDS, T_ITAG>
		>
{
public:
	template<typename TR_GEOM_TRAITS, typename TR_TDS, typename TR_ITAG>
	struct Rebind {
		using type = LIB_DTS2_NAMESPACE::Constrained_delaunay_triangulation_s2<TR_GEOM_TRAITS, TR_TDS, TR_ITAG>;
	};
public:
	using MyBaseClass = Constrained_triangulation_base_s2< CGAL::Constrained_Delaunay_triangulation_2<T_GEOM_TRAITS, T_TDS, T_ITAG> >;
public:
	using Triangulation_data_structure = typename MyBaseClass::Triangulation_data_structure;
	using Intersection_tag = typename MyBaseClass::Intersection_tag;
	using Geom_traits = typename MyBaseClass::Geom_traits;
	using Trs = typename MyBaseClass::Trs;
	using Tds = typename MyBaseClass::Tds;
	using Itag = typename MyBaseClass::Itag;
public:
	using FT = typename MyBaseClass::FT;
	using Point = typename MyBaseClass::Point;
	using Point_3 = typename MyBaseClass::Point_3;
	using GeoCoord = typename MyBaseClass::GeoCoord;
	using SphericalCoord = typename MyBaseClass::SphericalCoord;
public: //triangulation types
	using Segment = typename MyBaseClass::Segment;
	using Edge = typename MyBaseClass::Edge;
	
	using Vertex = typename MyBaseClass::Vertex;
	using Vertex_handle = typename MyBaseClass::Vertex_handle;
	using Vertex_circulator = typename MyBaseClass::Vertex_circulator;

	using Face = typename MyBaseClass::Face;
	using Face_handle = typename MyBaseClass::Face_handle;
	using Face_circulator = typename MyBaseClass::Face_circulator;
	
	using Edge_circulator = typename MyBaseClass::Edge_circulator;
	
	using Finite_vertices_iterator = typename MyBaseClass::Finite_vertices_iterator;
	using Finite_edges_iterator = typename MyBaseClass::Finite_edges_iterator;
	using Finite_faces_iterator = typename MyBaseClass::Finite_faces_iterator;
public:
	using Project_on_sphere = typename MyBaseClass::Project_on_sphere;

public:
	Constrained_delaunay_triangulation_s2(const Geom_traits & traits);
	Constrained_delaunay_triangulation_s2(int precision);
	Constrained_delaunay_triangulation_s2(Constrained_delaunay_triangulation_s2 && other);
	Constrained_delaunay_triangulation_s2(const Constrained_delaunay_triangulation_s2 & other) = delete;
	Constrained_delaunay_triangulation_s2 & operator=(Constrained_delaunay_triangulation_s2 && other);
	
public: //insertion
	using MyBaseClass::insert;
	void insert(const Vertex_handle & a, const Vertex_handle & b);
	void insert(const Point & a, const Point & b, bool snap = false);
	void insert(const GeoCoord & a, const GeoCoord & b);
	void insert(const SphericalCoord & a, const SphericalCoord & b);

	template<typename T_POINT_ITERATOR, typename T_INDICES_ITERATOR>
	void insert_constraints(T_POINT_ITERATOR points_first, T_POINT_ITERATOR points_last, T_INDICES_ITERATOR indices_first, T_INDICES_ITERATOR indices_last);
public: //info
	bool is_constrained(const Edge & e) const;
};

///The following typedefs may be of interest for users

template<typename T_TDS =
	typename internal::ConstrainedTriangulationDataStructureSelector<
		Constrained_delaunay_triangulation_traits_s2<CGAL::Exact_predicates_exact_constructions_kernel>
	>
>
using Constrained_Delaunay_triangulation_no_intersections_s2 =
	Constrained_delaunay_triangulation_s2<
		typename T_TDS::Geom_traits,
		T_TDS, CGAL::No_intersection_tag>;

template<typename T_TDS = typename internal::ConstrainedTriangulationDataStructureSelector<Constrained_delaunay_triangulation_with_inexact_intersections_traits_s2> >
using Constrained_Delaunay_triangulation_with_inexact_intersections_s2 =
	Constrained_delaunay_triangulation_s2<
		Constrained_delaunay_triangulation_with_inexact_intersections_traits_s2,
		T_TDS, CGAL::Exact_predicates_tag>;

template<typename T_TDS = typename internal::ConstrainedTriangulationDataStructureSelector<Constrained_delaunay_triangulation_with_inexact_intersections_traits_s2_64> >
using Constrained_Delaunay_triangulation_with_inexact_intersections_s2_64 =
	Constrained_delaunay_triangulation_s2<
		Constrained_delaunay_triangulation_with_inexact_intersections_traits_s2_64,
		T_TDS, CGAL::Exact_predicates_tag>;
		
template<typename T_TDS = typename internal::ConstrainedTriangulationDataStructureSelector<Constrained_delaunay_triangulation_with_inexact_intersections_traits_s2_sp> >
using Constrained_Delaunay_triangulation_with_inexact_intersections_s2_sp =
	Constrained_delaunay_triangulation_s2<
		Constrained_delaunay_triangulation_with_inexact_intersections_traits_s2_sp,
		T_TDS, CGAL::Exact_predicates_tag>;

template<typename T_TDS = typename internal::ConstrainedTriangulationDataStructureSelector<Constrained_delaunay_triangulation_with_exact_intersections_traits_s2> >
using Constrained_Delaunay_triangulation_with_exact_intersections_s2 =
	Constrained_delaunay_triangulation_s2<
		Constrained_delaunay_triangulation_with_exact_intersections_traits_s2,
		T_TDS, CGAL::Exact_intersections_tag>;
		
template<typename T_TDS = typename internal::ConstrainedTriangulationDataStructureSelector<Constrained_delaunay_triangulation_with_exact_intersections_spherical_traits_s2> >
using Constrained_Delaunay_triangulation_with_exact_intersections_spherical_s2 =
	Constrained_delaunay_triangulation_s2<
		Constrained_delaunay_triangulation_with_exact_intersections_spherical_traits_s2,
		T_TDS, CGAL::Exact_intersections_tag>;

template<typename T_VERTEX_INFO, typename T_FACE_INFO>
using Constrained_Delaunay_triangulation_no_intersections_with_info_s2 =
	Constrained_delaunay_triangulation_s2<
		Constrained_delaunay_triangulation_traits_s2<CGAL::Exact_predicates_exact_constructions_kernel>,
		typename internal::ConstrainedTriangulationDataStructureSelector<
			Constrained_delaunay_triangulation_traits_s2<CGAL::Exact_predicates_exact_constructions_kernel>,
			T_VERTEX_INFO,
			T_FACE_INFO
		>::type,
		CGAL::No_intersection_tag
	>;

template<typename T_VERTEX_INFO, typename T_FACE_INFO>
using Constrained_Delaunay_triangulation_with_inexact_intersections_with_info_s2 =
	Constrained_delaunay_triangulation_s2<
		Constrained_delaunay_triangulation_with_inexact_intersections_traits_s2,
		typename internal::ConstrainedTriangulationDataStructureSelector<
			Constrained_delaunay_triangulation_with_inexact_intersections_traits_s2,
			T_VERTEX_INFO,
			T_FACE_INFO
		>::type,
		CGAL::Exact_predicates_tag
	>;

template<typename T_VERTEX_INFO, typename T_FACE_INFO>
using Constrained_Delaunay_triangulation_with_inexact_intersections_with_info_s2_64 =
	Constrained_delaunay_triangulation_s2<
		Constrained_delaunay_triangulation_with_inexact_intersections_traits_s2_64,
		typename internal::ConstrainedTriangulationDataStructureSelector<
			Constrained_delaunay_triangulation_with_inexact_intersections_traits_s2_64,
			T_VERTEX_INFO,
			T_FACE_INFO
		>::type,
		CGAL::Exact_predicates_tag
	>;
	
template<typename T_VERTEX_INFO, typename T_FACE_INFO>
using Constrained_Delaunay_triangulation_with_inexact_intersections_with_info_s2_sp =
	Constrained_delaunay_triangulation_s2<
		Constrained_delaunay_triangulation_with_inexact_intersections_traits_s2_sp,
		typename internal::ConstrainedTriangulationDataStructureSelector<
			Constrained_delaunay_triangulation_with_inexact_intersections_traits_s2_sp,
			T_VERTEX_INFO,
			T_FACE_INFO
		>::type,
		CGAL::Exact_predicates_tag
	>;

template<typename T_VERTEX_INFO, typename T_FACE_INFO>
using Constrained_Delaunay_triangulation_with_inexact_intersections_with_info_s2_spk =
	Constrained_delaunay_triangulation_s2<
		Kernel_sp,
		typename internal::ConstrainedTriangulationDataStructureSelector<
			Kernel_sp,
			T_VERTEX_INFO,
			T_FACE_INFO
		>::type,
		CGAL::Exact_predicates_tag
	>;
	
template<typename T_VERTEX_INFO, typename T_FACE_INFO>
using Constrained_Delaunay_triangulation_with_inexact_intersections_with_info_s2_spk_64 =
	Constrained_delaunay_triangulation_s2<
		Kernel_sp_64,
		typename internal::ConstrainedTriangulationDataStructureSelector<
			Kernel_sp_64,
			T_VERTEX_INFO,
			T_FACE_INFO
		>::type,
		CGAL::Exact_predicates_tag
	>;

template<typename T_VERTEX_INFO, typename T_FACE_INFO>
using Constrained_Delaunay_triangulation_with_exact_intersections_with_info_s2 =
	Constrained_delaunay_triangulation_s2<
		Constrained_delaunay_triangulation_with_exact_intersections_traits_s2,
		typename internal::ConstrainedTriangulationDataStructureSelector<
			Constrained_delaunay_triangulation_with_exact_intersections_traits_s2,
			T_VERTEX_INFO,
			T_FACE_INFO
		>::type,
		CGAL::Exact_intersections_tag
	>;
	
template<typename T_VERTEX_INFO, typename T_FACE_INFO>
using Constrained_Delaunay_triangulation_with_exact_intersections_spherical_with_info_s2 =
	Constrained_delaunay_triangulation_s2<
		Constrained_delaunay_triangulation_with_exact_intersections_spherical_traits_s2,
		typename internal::ConstrainedTriangulationDataStructureSelector<
			Constrained_delaunay_triangulation_with_exact_intersections_spherical_traits_s2,
			T_VERTEX_INFO,
			T_FACE_INFO
		>::type,
		CGAL::Exact_intersections_tag
	>;

}//end namespace LIB_DTS2_NAMESPACE

//now the implementation
namespace LIB_DTS2_NAMESPACE {

#define TMPL_HDR template<typename T_GEOM_TRAITS, typename T_TDS, typename T_ITAG>
#define TMPL_CLS Constrained_delaunay_triangulation_s2<T_GEOM_TRAITS, T_TDS, T_ITAG>

TMPL_HDR
TMPL_CLS::Constrained_delaunay_triangulation_s2(const Geom_traits & traits) :
MyBaseClass(traits)
{}
	
TMPL_HDR
TMPL_CLS::Constrained_delaunay_triangulation_s2(int _precision) :
MyBaseClass(_precision)
{}

TMPL_HDR
TMPL_CLS::Constrained_delaunay_triangulation_s2(TMPL_CLS && other) :
MyBaseClass( std::move(other) )
{}

TMPL_HDR
TMPL_CLS &
TMPL_CLS::operator=(TMPL_CLS && other) {
	MyBaseClass::operator=( std::move(other) );
	return *this;
}

TMPL_HDR
bool
TMPL_CLS::is_constrained(const Edge& e) const {
	return MyBaseClass::trs().is_constrained(e);
}

TMPL_HDR
void
TMPL_CLS::insert(const Vertex_handle & a, const Vertex_handle & b) {
	MyBaseClass::trs().insert(a, b);
}

TMPL_HDR
void
TMPL_CLS::insert(const Point & a, const Point & b, bool snap) {
	if (snap) {
		Point_3 sa = MyBaseClass::proj()(a);
		Point_3 sb = MyBaseClass::proj()(b);
		MyBaseClass::trs().insert(sa, sb);
	}
	else {
		MyBaseClass::trs().insert(a, b);
	}
}

TMPL_HDR
void
TMPL_CLS::insert(const GeoCoord & a, const GeoCoord & b) {
	Point_3 sa = MyBaseClass::proj()(a);
	Point_3 sb = MyBaseClass::proj()(b);
	MyBaseClass::trs().insert(sa, sb);
}

TMPL_HDR
void
TMPL_CLS::insert(const SphericalCoord & a, const SphericalCoord & b) {
	Point_3 sa = MyBaseClass::proj()(a);
	Point_3 sb = MyBaseClass::proj()(b);
	MyBaseClass::trs().insert(sa, sb);
}

TMPL_HDR
template<typename T_POINT_ITERATOR, typename T_INDICES_ITERATOR>
void
TMPL_CLS::insert_constraints(T_POINT_ITERATOR points_first, T_POINT_ITERATOR points_last, T_INDICES_ITERATOR indices_first, T_INDICES_ITERATOR indices_last) {
	std::vector<Point_3> tmp;
	using std::distance;
	tmp.reserve(distance(points_first, points_last));
	for(auto it(points_first); it != points_last; ++it) {
		tmp.emplace_back(MyBaseClass::proj()(*it));
	}
	MyBaseClass::trs().insert_constraints(tmp.begin(), tmp.end(), indices_first, indices_last);
}

#undef TMPL_CLS
#undef TMPL_HDR

}//end namespace LIB_DTS2_NAMESPACE

#endif
