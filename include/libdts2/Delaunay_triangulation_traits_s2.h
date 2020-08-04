#pragma once
#ifndef LIB_DTS2_DELAUNAY_TRIANGULATION_TRAITS_S2_H
#define LIB_DTS2_DELAUNAY_TRIANGULATION_TRAITS_S2_H

#include <libdts2/constants.h>
#include <libdts2/Constrained_delaunay_triangulation_traits_s2.h>

namespace LIB_DTS2_NAMESPACE {

template<typename T_LINEAR_KERNEL = CGAL::Exact_predicates_exact_constructions_kernel>
using Delaunay_triangulation_traits_s2 = Constrained_delaunay_triangulation_traits_s2<T_LINEAR_KERNEL>;

}//end namespace LIB_DTS2_NAMESPACE


#endif
