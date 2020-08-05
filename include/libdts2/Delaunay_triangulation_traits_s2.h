#pragma once
#ifndef LIB_DTS2_DELAUNAY_TRIANGULATION_TRAITS_S2_H
#define LIB_DTS2_DELAUNAY_TRIANGULATION_TRAITS_S2_H

#include <libdts2/constants.h>
#include <libdts2/Constrained_delaunay_triangulation_traits_s2.h>

namespace LIB_DTS2_NAMESPACE {

template<
	typename T_LINEAR_KERNEL,
	typename T_AUX_POINT_GENERATOR=detail::EpsBasedAuxPoints<T_LINEAR_KERNEL>,
	typename T_POINT = typename T_LINEAR_KERNEL::Point_3
>
using Delaunay_triangulation_traits_s2 = Constrained_delaunay_triangulation_traits_s2<T_LINEAR_KERNEL, T_AUX_POINT_GENERATOR, T_POINT>;

}//end namespace LIB_DTS2_NAMESPACE


#endif
