#pragma once
#ifndef LIB_DTS2_DEBUG_H
#define LIB_DTS2_DEBUG_H

#include <libdts2/Constrained_delaunay_triangulation_s2.h>

extern "C" {
	
void dts2_debug_print_point_sp(LIB_DTS2_NAMESPACE::Constrained_delaunay_triangulation_with_inexact_intersections_traits_s2_sp::Point_2 const & p);
void dts2_debug_print_kernel_sp_point(LIB_DTS2_NAMESPACE::Kernel_sp::Point_2 const & p);
void dts2_debug_print_kernel_sp64_point(LIB_DTS2_NAMESPACE::Kernel_sp_64::Point_2 const & p);

}

namespace LIB_DTS2_NAMESPACE {
	void init_interactive_debuging();
}

#endif
