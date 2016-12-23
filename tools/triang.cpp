
/**
  * Graph format is as follows
  * number_of_points
  * number_of_constraints
  * points
  * constraints
  *
  * Points may be encoded in any format accepted my ratss::InputPoint
  * Constraints reference points: id1 id2
  * 
  */

#include <fstream>
#include <libdts2/Constrained_delaunay_triangulation_s2.h>
#include <libratss/mpreal.h>
#include <libratss/GeoCoord.h>
#include <libratss/SphericalCoord.h>
#include "../vendor/libratss/common/points.h"

struct FloatPoint3 {
	mpfr::mpreal coords[3];
	mpfr::mpreal & x() { return coords[0]; }
	mpfr::mpreal & y() { return coords[1]; }
	mpfr::mpreal & z() { return coords[2]; }
};

typedef enum {TT_DELAUNAY, TT_CONSTRAINED, TT_CONSTRAINED_INEXACT, TT_CONSTRAINED_EXACT} TriangulationType;

struct VertexInfo {
	int id;
};

///vertices need to have VertexInfo as info
template<typename TRS>
struct TriangulationWriter {
	using Finite_vertices_iterator = typename TRS::Finite_vertices_iterator; 
	using Finite_edges_iterator = typename TRS::Finite_edges_iterator;
	using Point = typename TRS::Point;
	using FT = typename TRS::FT;
	
	ratss::OutputPoint::Format outFormat;
	
	void write(std::ostream & out, TRS & trs) {
		using std::distance;
		std::size_t vertexCount = distance(trs.finite_vertices_begin(), trs.finite_vertices_end());
		out << vertexCount << '\n';
		out << distance(trs.finite_edges_begin(), trs.finite_edges_end()) << '\n';
		
		ratss::OutputPoint op(3);
		std::vector<int> vertexId2OutId(vertexCount);
		std::size_t counter = 0;
		for(Finite_vertices_iterator it(trs.finite_vertices_begin()), end(trs.finite_vertices_end()); it != end; ++it) {
			const VertexInfo & vi = it->info();
			const Point & p = it->point();
			op.coords[0] = ratss::Conversion<FT>::toMpq( p.x() );
			op.coords[1] = ratss::Conversion<FT>::toMpq( p.y() );
			op.coords[2] = ratss::Conversion<FT>::toMpq( p.z() );
			
			op.print(out, outFormat);
			vertexId2OutId.at(vi.id) = counter;
			++counter;
		}
		
		for(Finite_edges_iterator it(trs.finite_edges_begin()), end(trs.finite_edges_end()) ; it != end; ++it) {
			auto & face = *(it->first);
			int index = it->second;
			auto vs = face.vertex(face.cw(index));
			auto vt = face.vertex(face.ccw(index));
			const VertexInfo & vis = vs->info();
			const VertexInfo & vit = vt->info();
			out << vis.id << ' ' << vit.id << '\n';
		}
	}
};

struct TriangulationCreator {
	virtual void create(const std::vector<FloatPoint3> & points, const std::vector<std::pair<int, int>> & edges) = 0;
	virtual void write(std::ostream & out, ratss::OutputPoint::Format format) = 0;
};

struct TriangulationCreatorDelaunay {
	using Tr =  dts2::Delaunay_triangulation_with_info_s2<VertexInfo, void>;
	
	virtual void create(const std::vector<FloatPoint3> & points, const std::vector<std::pair<int, int>> & edges) {
		
	}
	
	virtual void write(std::ostream & out, ratss::OutputPoint::Format format) {
		
	}
	
};

void help(std::ostream & out) {
	out << "triang OPTIONS:\n"
		"-t type\ttype = [d,delaunay, c,constrained,cx,constrained-intersection,cxe,constrained-intesection-exact]\n"
		"-p precision\tprecision to calculate intermediate results\n"
		"-e signficands\tprecision with which the points on the sphere are calculated\n"
		"-if format\tformat of input points\n"
		"-of format\tformat of output points\n"
	<< std::endl;
}

struct Config {
	TriangulationType triangType;
	ratss::InputPoint::Format inFormat;
	ratss::OutputPoint::Format outFormat;
	int precision;
	int signficands;
};

int main(int argc, char ** argv) {
	
}