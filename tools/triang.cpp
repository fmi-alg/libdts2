
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
#include <libratss/ProjectSN.h>
#include <libratss/util/BasicCmdLineOptions.h>
#include <libratss/util/InputOutputPoints.h>
#include <libratss/Conversion.h>

using K = CGAL::Exact_predicates_exact_constructions_kernel;
using Point3 = K::Point_3;

typedef enum {TT_DELAUNAY, TT_CONSTRAINED, TT_CONSTRAINED_INEXACT, TT_CONSTRAINED_EXACT} TriangulationType;

struct VertexInfo {
	VertexInfo() : id(-1) {}
	VertexInfo(int id) : id(id) {}
	bool valid() const { return id >= 0; }
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
	
	TriangulationWriter(ratss::OutputPoint::Format outFormat) : outFormat(outFormat) {}
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

class TriangulationCreator {
public:
	TriangulationCreator() {}
	virtual ~TriangulationCreator() {}
public:
	virtual void create(const std::vector<std::pair<Point3, VertexInfo>> & points, const std::vector<std::pair<int, int>> & edges) = 0;
	virtual void write(std::ostream & out, ratss::OutputPoint::Format format) = 0;
};

class TriangulationCreatorDelaunay: public TriangulationCreator {
public:
	using Tr =  dts2::Delaunay_triangulation_with_info_s2<VertexInfo, void>;
	using Point_3 = Tr::Point_3;
private:
	Tr m_tr;
public:
	TriangulationCreatorDelaunay(int significands) : m_tr(significands) {}
	
	virtual void create(const std::vector<std::pair<Point3, VertexInfo>> & points, const std::vector<std::pair<int, int>> & edges) override {
		m_tr.insert(points.begin(), points.end());
	}
	
	virtual void write(std::ostream & out, ratss::OutputPoint::Format format) override {
		TriangulationWriter<Tr> writer(format);
		writer.write(out, m_tr);
	}
};

template< template<typename, typename> class T_TRS>
class TriangulationCreatorConstrainedDelaunay: public TriangulationCreator {
public:
	using Tr =  T_TRS<VertexInfo, void>;
	using Vertex_handle = typename Tr::Vertex_handle;
	using Finite_vertices_iterator = typename Tr::Finite_vertices_iterator;
private:
	Tr m_tr;
public:
	TriangulationCreatorConstrainedDelaunay(int significands) : m_tr(significands) {}
	
	virtual void create(const std::vector<std::pair<Point3, VertexInfo>> & points, const std::vector<std::pair<int, int>> & edges) override  {
		m_tr.insert(points.begin(), points.end());
		std::vector<Vertex_handle> pId2Vertex(points.size());
		for(Finite_vertices_iterator it(m_tr.finite_vertices_begin()), end(m_tr.finite_vertices_end()); it != end; ++it) {
			const VertexInfo & vi = it->info();
			if (vi.valid()) {
				Vertex_handle & vh = pId2Vertex.at(vi.id);
				vh =(Vertex_handle) it;
			}
		}
		for(const std::pair<int, int> & e : edges) {
			m_tr.insert(pId2Vertex.at(e.first), pId2Vertex.at(e.second));
		}
	}
	
	virtual void write(std::ostream & out, ratss::OutputPoint::Format format) override  {
		TriangulationWriter<Tr> writer(format);
		writer.write(out, m_tr);
	}
};

using TriangulationCreatorNoIntersectionsConstrainedDelaunay = TriangulationCreatorConstrainedDelaunay<dts2::Constrained_Delaunay_triangulation_no_intersections_with_info_s2>;
using TriangulationCreatorInExactIntersectionsConstrainedDelaunay = TriangulationCreatorConstrainedDelaunay<dts2::Constrained_Delaunay_triangulation_with_inexact_intersections_with_info_s2>;

class TriangulationCreatorExactIntersectionsConstrainedDelaunay: public TriangulationCreator {
public:
	using Tr =  dts2::Constrained_Delaunay_triangulation_with_exact_intersections_with_info_s2<VertexInfo, void>;
	using Point_3 = Tr::Point_3;
private:
	Tr m_tr;
public:
	TriangulationCreatorExactIntersectionsConstrainedDelaunay(int significands) : m_tr(significands) {}
	
	virtual void create(const std::vector<std::pair<Point3, VertexInfo>> & points, const std::vector<std::pair<int, int>> & edges) override {
		std::vector< std::pair<Point_3, VertexInfo> > myPoints;
		myPoints.reserve(points.size());
		for(const std::pair<Point3, VertexInfo> & pi : points) {
			CORE::Expr x = ratss::Conversion<CORE::Expr>::moveFrom( ratss::Conversion<K::FT>::toMpq(pi.first.x()) );
			CORE::Expr y = ratss::Conversion<CORE::Expr>::moveFrom( ratss::Conversion<K::FT>::toMpq(pi.first.y()) );
			CORE::Expr z = ratss::Conversion<CORE::Expr>::moveFrom( ratss::Conversion<K::FT>::toMpq(pi.first.z()) );
			Point_3 myPoint(x, y, z);
			myPoints.emplace_back(std::move(myPoint), pi.second);
		}
		m_tr.insert(myPoints.begin(), myPoints.end());
	}
	
	virtual void write(std::ostream & out, ratss::OutputPoint::Format format) override {
		TriangulationWriter<Tr> writer(format);
		writer.write(out, m_tr);
	}
};

class Config: public ratss::BasicCmdLineOptions {
public:
	TriangulationType triangType;
public:
	using ratss::BasicCmdLineOptions::parse;
public:
	virtual bool parse(const std::string & currentToken,int & i, int argc, char ** argv) override;
	void help(std::ostream & out) const;
	void print(std::ostream & out) const;
};

struct Data {
	std::vector<std::pair<Point3, VertexInfo>> points;
	std::vector<std::pair<int, int>> edges;
	TriangulationCreator * tc;
	
	void init(TriangulationType tt, int significands);
	void read(std::istream& is, const Config& cfg);
	void create();
	void write(std::ostream & out, ratss::OutputPoint::Format outFormat);
};

class InputOutput {
public:
	std::istream & input();
	std::ostream & output();
	void setInput(const std::string & inFileName);
	void setOutput(const std::string & outFileName);
private:
	std::istream * inFile = 0;
	std::ostream * outFile = 0;
	
	std::ifstream inFileHandle;
	std::ofstream outFileHandle;
	
};

///now the main

int main(int argc, char ** argv) {
	Config cfg;
	InputOutput io;
	Data data;
	
	int parseResult = cfg.parse(argc, argv);
	if(parseResult <= 0) {
		cfg.help(std::cerr);
		return parseResult;
	}
	
	cfg.print(std::cerr);
	
	io.setInput(cfg.inFileName);
	io.setOutput(cfg.outFileName);
	
	
	data.init(cfg.triangType, cfg.significands);
	data.read(io.input(), cfg);
	data.create();
	data.write(io.output(), cfg.outFormat);
	
}

//now the implementations of the functions of the classes


bool Config::parse(const std::string & token,int & i, int argc, char ** argv) {
	if (token == "-t" && i+1 < argc) {
		std::string type( argv[i+1] );
		if (type == "d" || type == "delaunay") {
			triangType = TT_DELAUNAY;
		}
		else if (type == "c" || type == "constrained") {
			triangType = TT_CONSTRAINED;
		}
		else if (type == "cx" || type == "constrained-intersection") {
			triangType = TT_CONSTRAINED_INEXACT;
		}
		else if (type == "cxe" || type == "constrained-intersection-exact") {
			triangType = TT_CONSTRAINED_EXACT;
		}
		else {
			throw ratss::BasicCmdLineOptions::ParseError("Unknown triangulation type: " + type);
		}
		++i;
	}
	else {
		return false;
	}
	return true;
}

void Config::help(std::ostream & out) const {
	out << "triang OPTIONS:\n"
		"-t type\ttype = [d,delaunay, c,constrained,cx,constrained-intersection,cxe,constrained-intesection-exact]\n";
	ratss::BasicCmdLineOptions::options_help(out);
	out << std::endl;
}


void Config::print(std::ostream & out) const {
	out << "Triangulation type: ";
	switch (triangType) {
	case TT_DELAUNAY:
		out << "delaunay";
		break;
	case TT_CONSTRAINED:
		out << "constrained no intersections";
		break;
	case TT_CONSTRAINED_INEXACT:
		out << "constrained in-exact intersecionts";
		break;
	case TT_CONSTRAINED_EXACT:
		out << "constrained exact intersections";
		break;
	default:
		out << "invalid";
		break;
	}
	out << '\n';
	ratss::BasicCmdLineOptions::options_selection(out);
}

std::istream & InputOutput::input() {
	return *inFile;
}

std::ostream & InputOutput::output() {
	return *outFile;
}

void InputOutput::setInput(const std::string & inFileName) {
	assert(!inFile);
	if (inFileName.size()) {
		inFileHandle.open(inFileName);
		if (!inFileHandle.is_open()) {
			throw std::runtime_error("Could not open input file: " + inFileName);
		}
		inFile = &inFileHandle;
	}
	else {
		inFile = &std::cin;
	}

}

void InputOutput::setOutput(const std::string & outFileName) {
	assert(!outFile);
	if (outFileName.size()) {
		outFileHandle.open(outFileName);
		if (!outFileHandle.is_open()) {
			throw std::runtime_error("Could not open output file: " + outFileName);
		}
		outFile = &outFileHandle;
	}
	else {
		outFile = &std::cout;
	}
}

void Data::init(TriangulationType tt, int significands) {
	TriangulationCreator * tc = 0;
	switch (tt) {
	case TT_DELAUNAY:
		tc = new TriangulationCreatorDelaunay(significands);
		break;
	case TT_CONSTRAINED:
		tc = new TriangulationCreatorNoIntersectionsConstrainedDelaunay(significands);
		break;
	case TT_CONSTRAINED_INEXACT:
		tc = new TriangulationCreatorInExactIntersectionsConstrainedDelaunay(significands);
		break;
	case TT_CONSTRAINED_EXACT:
		tc = new TriangulationCreatorExactIntersectionsConstrainedDelaunay(significands);
		break;
	default:
		throw std::runtime_error("Unkown triangulation type");
	};
}

void Data::read(std::istream& is, const Config & cfg) {
	std::size_t num_points, num_edges;
	is >> num_points >> num_edges;
	
	std::cerr << "Need to fetch " << num_points << " points and " << num_edges << " edges" << std::endl;
	
	points.reserve(num_points);
	edges.reserve(num_edges);
	
	ratss::InputPoint ip;
	ratss::OutputPoint op;
	ratss::ProjectSN proj;
	
	if (cfg.progress) {
		std::cerr << std::endl;
	}
	std::size_t counter = 0;
	while( is.good() && points.size() < num_points) {
		if (is.peek() == '\n') {
			is.get();
			continue;
		}
		ip.assign(is, cfg.inFormat, cfg.precision);
		if (cfg.normalize) {
			ip.normalize();
		}
		ip.setPrecision(cfg.precision);
		op.clear();
		op.resize(ip.coords.size());
		proj.snap(ip.coords.begin(), ip.coords.end(), op.coords.begin(), cfg.snapType, cfg.significands);
		
		K::FT x = ratss::Conversion<K::FT>::moveFrom(op.coords.at(0));
		K::FT y = ratss::Conversion<K::FT>::moveFrom(op.coords.at(1));
		K::FT z = ratss::Conversion<K::FT>::moveFrom(op.coords.at(2));
		Point3 p3(x, y, z);
		VertexInfo vi(counter);
		points.emplace_back(std::move(p3), vi);

		++counter;
		if (cfg.progress && counter % 1000 == 0) {
			std::cerr << '\xd' << counter/1000 << "k" << std::flush;
		}
	}
	//now read in the edges
	while( is.good() && edges.size() < num_edges) {
		if (is.peek() == '\n') {
			is.get();
			continue;
		}
		std::pair<int, int> e;
		is >> e.first >> e.second;
		edges.emplace_back(e);
	}
	
	if (num_points != points.size()) {
		throw std::runtime_error("Read " + std::to_string(points.size()) + " points instead of " + std::to_string(num_points));
	}
	if (num_edges != edges.size()) {
		throw std::runtime_error("Read " + std::to_string(edges.size()) + " edges instead of " + std::to_string(num_edges));
	}
}

void Data::create() {
	tc->create(points, edges);
	points.clear();
	edges.clear();
}

void Data::write(std::ostream & out, ratss::OutputPoint::Format outFormat) {
	tc->write(out, outFormat);
}