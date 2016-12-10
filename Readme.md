This is a small library to compute a delaunay triangulation of points that are on the sphere.
If input points are not exactly on the sphere then they are snaped on the sphere using libratss.
It essentialy provides an adaptor class for the (Delaunay) triangulation algorithms that are provides by CGAL.
Any other triangulation algorithms/data structures following the triangulation concepts of CGAL should work aswell.

This library currently supports Delaunay and constrained Delaunay-triangulations without intersecting segments.
There is untested code to also support constrained Delaunay-triangulations where segments are allowed to intersect.
There is the possibility to use exact arithmetic to compute the intersecting point using CORE::Expr.
Approximate intersections are of course faster.
There is also preliminary support tu use the spherical kernel as a base.