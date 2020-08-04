#pragma once
#ifndef LIB_DTS2_DELAUNAY_TRIANGULATION_S2_H
#define LIB_DTS2_DELAUNAY_TRIANGULATION_S2_H

#include <libdts2/constants.h>
#include <libdts2/Triangulation_base_s2.h>

//we have to include these first (due to some overloaded internal functions)
#include <libdts2/Delaunay_triangulation_traits_s2.h>

#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Constrained_triangulation_face_base_2.h>
#include <CGAL/Triangulation_data_structure_2.h>

#include <CGAL/Delaunay_triangulation_2.h>

namespace LIB_DTS2_NAMESPACE {
namespace internal {
	template<typename T_KERNEL, typename T_INFO = void>
	struct VertexBaseSelector {
		using type = CGAL::Triangulation_vertex_base_with_info_2<T_INFO, T_KERNEL>;
	};

	template<typename T_KERNEL, typename T_INFO = void>
	struct FaceBaseSelector {
		using type = CGAL::Triangulation_face_base_with_info_2<T_INFO, T_KERNEL>;
	};

	template<typename T_KERNEL>
	struct VertexBaseSelector<T_KERNEL, void> {
		using type = CGAL::Triangulation_vertex_base_2<T_KERNEL>;
	};

	template<typename T_KERNEL>
	struct FaceBaseSelector<T_KERNEL, void> {
		using type = CGAL::Triangulation_face_base_2<T_KERNEL>;
	};
	
	template<typename T_KERNEL, typename T_VERTEX_INFO=void, typename T_FACE_INFO=void>
	class TriangulationDataStructureSelector {
	private:
		using Vb = typename VertexBaseSelector<T_KERNEL, T_VERTEX_INFO>::type;
		using Fb = typename FaceBaseSelector<T_KERNEL, T_FACE_INFO>::type;
	public:
		using type = CGAL::Triangulation_data_structure_2<Vb, Fb>;
	};
} //end namespace internal


/** This is a delaunay triangulation living on the sphere S2
	It uses CGAL::Delaunay_triangulation_2 as bases
**/
template<typename T_GEOM_TRAITS, typename T_TDS>
class Delaunay_triangulation_s2:
	public Triangulation_base_s2<
		CGAL::Delaunay_triangulation_2<T_GEOM_TRAITS, T_TDS>
		>
{
public:
	using MyBaseClass = Triangulation_base_s2< CGAL::Delaunay_triangulation_2<T_GEOM_TRAITS, T_TDS> >;
public:
	using Triangulation_data_structure = typename MyBaseClass::Triangulation_data_structure;
	using Geom_traits = typename MyBaseClass::Geom_traits;
	using Trs = typename MyBaseClass::Trs;
	using Tds = typename MyBaseClass::Tds;
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
	Delaunay_triangulation_s2(Geom_traits const & traits);
	Delaunay_triangulation_s2(int significands);
	Delaunay_triangulation_s2(Delaunay_triangulation_s2 && other);
	Delaunay_triangulation_s2(const Delaunay_triangulation_s2 & other) = delete;
	Delaunay_triangulation_s2 & operator=(Delaunay_triangulation_s2 && other);
};


template<typename T_VERTEX_INFO, typename T_FACE_INFO, typename T_KERNEL = CGAL::Exact_predicates_exact_constructions_kernel>
using Delaunay_triangulation_with_info_s2 =
	Delaunay_triangulation_s2<
		Delaunay_triangulation_traits_s2<T_KERNEL>,
		typename internal::TriangulationDataStructureSelector<
			Delaunay_triangulation_traits_s2<T_KERNEL>,
			T_VERTEX_INFO,
			T_FACE_INFO
		>::type
	>;

}//end namespace LIB_DTS2_NAMESPACE

//now the implementation
namespace LIB_DTS2_NAMESPACE {

#define TMPL_HDR template<typename T_GEOM_TRAITS, typename T_TDS>
#define TMPL_CLS Delaunay_triangulation_s2<T_GEOM_TRAITS, T_TDS>

TMPL_HDR
TMPL_CLS::Delaunay_triangulation_s2(Geom_traits const & traits) :
MyBaseClass(traits)
{}

TMPL_HDR
TMPL_CLS::Delaunay_triangulation_s2(int _significands) :
MyBaseClass(_significands)
{}

TMPL_HDR
TMPL_CLS::Delaunay_triangulation_s2(TMPL_CLS && other) :
MyBaseClass( std::move(other) )
{}

TMPL_HDR
TMPL_CLS &
TMPL_CLS::operator=(TMPL_CLS && other) {
	MyBaseClass::operator=( std::move(other) );
	return *this;
}

#undef TMPL_CLS
#undef TMPL_HDR

}//end namespace LIB_DTS2_NAMESPACE

#endif
