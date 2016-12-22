
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
#include <vendor/libratss/common/points.h>

struct FloatPoint3 {
	mpfr::mpreal x;
	mpfr::mpreal y;
	mpfr::mpreal z;
};

typedef enum {TT_DELAUNAY, TT_CONSTRAINED, TT_CONSTRAINED_INEXACT, TT_CONSTRAINED_EXACT} TriangulationType;

struct TriangulationCreator {
	virtual void create(const std::vector<FloatPoint3> & points, const std::vector<std::pair<int, int>> & edges) = 0;
	virtual void write(std::ostream & out, ratss::OutputPoint::Format format);
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
	TriangulationCreator triangType;
	ratss::InputPoint::Format inFormat;
	ratss::OutputPoint::Format outFormat;
	int precision;
	int signficands;
};

int main(int argc, char ** argv) {
	
}