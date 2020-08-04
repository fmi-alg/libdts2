#pragma once
#ifndef LIB_DTS2_KERNEL_SP_BASE_H
#define LIB_DTS2_KERNEL_SP_BASE_H

#include <CGAL/Simple_cartesian.h>
#include <libdts2/Kernel_sp/Point_sp.h>

namespace CGAL {

template < typename K_, typename FT_>
struct Kernel_sp_base
{
    typedef K_                                          Kernel;
    typedef FT_                                         FT;
    typedef Cartesian_base<K_,FT_>                      Self;
    typedef Cartesian_tag                               Rep_tag;
    typedef Cartesian_tag                               Kernel_tag;

    enum { Has_filtered_predicates = false };
    typedef Boolean_tag<Has_filtered_predicates> Has_filtered_predicates_tag;

    typedef CGAL::Object                                Object_2;
    typedef CGAL::Object                                Object_3;

    // Boolean   had originally been Bool. It was renamed to avoid a conflict
    // between a macro defined in Xlib.h poorly chosen to have the same name,
    // that is 'Bool'.
    typedef typename Same_uncertainty_nt<bool, FT>::type
                                                        Boolean;
    typedef typename Same_uncertainty_nt<CGAL::Sign, FT>::type
                                                        Sign;
    typedef typename Same_uncertainty_nt<CGAL::Comparison_result, FT>::type
                                                        Comparison_result;
    typedef typename Same_uncertainty_nt<CGAL::Orientation, FT>::type
                                                        Orientation;
    typedef typename Same_uncertainty_nt<CGAL::Oriented_side, FT>::type
                                                        Oriented_side;
    typedef typename Same_uncertainty_nt<CGAL::Bounded_side, FT>::type
                                                        Bounded_side;
    typedef typename Same_uncertainty_nt<CGAL::Angle, FT>::type
                                                        Angle;

    template <typename T>
    struct Ambient_dimension {
      typedef typename T::Ambient_dimension type;
    };

    template <typename T>
    struct Feature_dimension {
      typedef typename T::Feature_dimension type;
    };

    typedef Point_sp<Kernel>                            Point_2;
    typedef VectorC2<Kernel>                            Vector_2;
    typedef DirectionC2<Kernel>                         Direction_2;
    typedef SegmentC2<Kernel>                           Segment_2;
    typedef LineC2<Kernel>                              Line_2;
    typedef RayC2<Kernel>                               Ray_2;
    typedef TriangleC2<Kernel>                          Triangle_2;
    typedef CircleC2<Kernel>                            Circle_2;
    typedef Iso_rectangleC2<Kernel>                     Iso_rectangle_2;
    typedef Aff_transformationC2<Kernel>                Aff_transformation_2;

    typedef PointC3<Kernel>                             Point_3;
    typedef VectorC3<Kernel>                            Vector_3;
    typedef DirectionC3<Kernel>                         Direction_3;
    typedef LineC3<Kernel>                              Line_3;
    typedef PlaneC3<Kernel>                             Plane_3;
    typedef RayC3<Kernel>                               Ray_3;
    typedef SegmentC3<Kernel>                           Segment_3;
    typedef TriangleC3<Kernel>                          Triangle_3;
    typedef TetrahedronC3<Kernel>                       Tetrahedron_3;
    typedef Iso_cuboidC3<Kernel>                        Iso_cuboid_3;
    typedef SphereC3<Kernel>                            Sphere_3;
    typedef CircleC3<Kernel>                            Circle_3;
    typedef Aff_transformationC3<Kernel>                Aff_transformation_3;

    typedef typename cpp11::array<FT_, 2>::const_iterator Cartesian_const_iterator_2;
    typedef typename cpp11::array<FT_, 3>::const_iterator Cartesian_const_iterator_3;

    // Undocumented stuff.
    typedef Data_accessorC2<Kernel>                     Data_accessor_2;
    typedef ConicCPA2<Point_2,Data_accessor_2>          Conic_2;
};

} //end namespace CGAL

#endif
