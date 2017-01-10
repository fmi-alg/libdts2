## Introduction
This is a small library to compute a delaunay triangulation of points that are on the sphere.
If input points are not exactly on the sphere then they are snaped on the sphere using libratss.
It essentialy provides an adaptor class for the (Delaunay) triangulation algorithms that are provided by CGAL.
Any other triangulation algorithms/data structures following the triangulation concepts of CGAL should work aswell.

## Scope
This library currently supports Delaunay and Constrained Delaunay Triangulations without intersecting segments.
There is also support for Constrained Delaunay Triangulations where segments are allowed to intersect.
There is the possibility to use exact arithmetic to compute the intersecting point using CORE::Expr.
Approximate intersections are of course faster.

### Things to consider
Computing points that are exactly on the sphere is hard. libdts2 relies on libratss to compute points on the sphere given a point that is near the sphere. Usually one has to deal with geographic coordinates or spherical coordinates. These have to be mapped to cartesian coordinates in 3D. Since floating-point computations are inexact this will not yield a point that is on the sphere. The distance to the sphere depends on the precision of the computation. The higher the precision, the closer the point to the sphere, the closer the snapped point to the original point.

Usually a computation precision of 53 is more than enough to represent points on a sphere representing our earth. This translates to a Point occupying 3 * 2 * 8 Bytes.

## Dependecies
* CGAL
* gmpxx

## Quickstart guide
### Building a triangulation
1. git clone --recursive https://github.com/fmi-alg/libdts2.git libdts2
2. cd libdts2 && mkdir build && cd build
3. cmake -DCMAKE_BUILD_TYPE=lto ../
4. cd dts2tools && make
5. ./triang -h

#### Getting sample graphs
There is a tool to create sample graphs out of OpenStreetMap data:

1. Download and build https://github.com/dbahrdt/OsmGraphCreator.git
2. Create graph of the type topotext and add the --no-reverse-edge switch
3. For triang select geo as an input type: -if geo

#### Visualizing graphs
There is a tool to visualize graphs create by the triang tool:

1. Download and build the dev-branch of https://github.com/invor/simplestGraphRendering
2. create an output of type -g simplest_andre
3. Start simplestGraphRendering with the options -f sg

### Use the library
1. add the repository to your project, make sure to checkout the subprojects of libdts2
2. make sure to add the appropriate compile definitions to your project. If you're using CMake and you're adding libdts2 as a subdirectory then these are exported. See the CMakeLists.txt for a list of exported variables.
3. add <libdts/Constrained_delaunay_triangulation_s2.h> as include
4. select an appropriate triangulation class, see below for a small example

#### Input formats
Currently input points can be any of:
* ratss::GeoCoord provides a small wrapper to handle wgs84 coordinates
* ratss::SphericalCoord provides a small wrapper to handle spherical coordinates
* Triangulation::Point provides a 3d Point

#### A small example

```C++
//This file has predefined typedefs for various triangulations.
#include <libdts/Constrained_delaunay_triangulation_s2.h>

//The following triangulation supports in-exact intersections and can store information in its vertices and faces:

template<typename VertexInfo = void, typename FaceInfo = void>
using MyTriangulation = Constrained_Delaunay_triangulation_with_inexact_intersections_with_info_s2<VertexInfo, FaceInfo>;

int significands = 53;
MyTriangulation tr(significands);
auto vh1 = tr.insert(ratss::GeoCoord(48.0, 8.0));
auto vh2 = tr.insert(ratss::GeoCoord(49.0, 9.0));
auto vh3 = tr.insert(ratss::GeoCoord(48.0, 9.0));
auto vh4 = tr.insert(ratss::GeoCoord(49.0, 8.0));
tr.insert(vh1, vh2);
tr.insert(vh3, vh4);
```

# License
libdts2 is licensed under LGPL v2.1
Since libdts2 links with mpfr, gmp and cgal the resulting binaries are usually governed by the GPL v2
