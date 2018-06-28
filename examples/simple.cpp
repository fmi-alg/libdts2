//This file has predefined typedefs for various triangulations.
#include <libdts2/Constrained_delaunay_triangulation_s2.h>

using namespace dts2;
//The following triangulation supports in-exact intersections and can store information in its vertices and faces:

template<typename VertexInfo = void, typename FaceInfo = void>
using MyTriangulation = Constrained_Delaunay_triangulation_with_inexact_intersections_with_info_s2<VertexInfo, FaceInfo>;


int main() {

	int significands = 53;
	MyTriangulation<void, void> tr(significands);
	auto vh1 = tr.insert(ratss::GeoCoord(48.0, 8.0));
	auto vh2 = tr.insert(ratss::GeoCoord(49.0, 9.0));
	auto vh3 = tr.insert(ratss::GeoCoord(48.0, 9.0));
	auto vh4 = tr.insert(ratss::GeoCoord(49.0, 8.0));
	tr.insert(vh1, vh2);
	tr.insert(vh3, vh4);
}
