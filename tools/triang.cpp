
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
#include <libratss/util/InputOutput.h>
#include <libratss/Conversion.h>

using K = CGAL::Exact_predicates_exact_constructions_kernel;
using Point3 = K::Point_3;

typedef enum {TT_DELAUNAY, TT_CONSTRAINED, TT_CONSTRAINED_INEXACT, TT_CONSTRAINED_EXACT, TT_CONSTRAINED_EXACT_SPHERICAL} TriangulationType;
typedef enum {GOT_INVALID, GOT_WITHOUT_SPECIAL, GOT_SIMPLEST_GRAPH_RENDERING, GOT_SIMPLEST_GRAPH_RENDERING_ANDRE} GraphOutputType;

struct VertexInfo {
	VertexInfo() : id(-1) {}
	VertexInfo(int id) : id(id) {}
	bool valid() const { return id >= 0; }
	int id;
};

bool is_constrained(const dts2::Delaunay_triangulation_with_info_s2<VertexInfo, void> & trs, const dts2::Delaunay_triangulation_with_info_s2<VertexInfo, void>::Edge & e) {
	return false;
}

template<typename T_TRS>
bool is_constrained(const T_TRS & trs, const typename T_TRS::Edge & e) {
	return trs.is_constrained(e);
}

using InputOutput = ratss::InputOutput;

///vertices need to have VertexInfo as info
template<typename TRS>
struct TriangulationWriter {
	using Vertex_handle = typename TRS::Vertex_handle;
	using Finite_vertices_iterator = typename TRS::Finite_vertices_iterator; 
	using Finite_edges_iterator = typename TRS::Finite_edges_iterator;
	using Finite_faces_iterator = typename TRS::Finite_faces_iterator;
	
	using Edge = typename TRS::Edge;

	using Point = typename TRS::Point;
	using FT = typename TRS::FT;
	
	ratss::OutputPoint::Format outFormat;
	GraphOutputType got;
	
	TriangulationWriter(ratss::OutputPoint::Format outFormat, GraphOutputType got) : outFormat(outFormat), got(got) {}
	
	void write(std::ostream & out, TRS & trs) {
		switch (got) {
		case GOT_WITHOUT_SPECIAL:
			writeWithoutSpecial(out, trs);
			break;
		case GOT_SIMPLEST_GRAPH_RENDERING:
			writeSimplestGraphRendering(out, trs);
			break;
		case GOT_SIMPLEST_GRAPH_RENDERING_ANDRE:
			writeSimplestGraphRenderingAndre(out, trs);
			break;
		default:
			throw std::runtime_error("Unknown graph output type: " + std::to_string(got));
		};
	}
	
	void writeWithoutSpecial(std::ostream & out, TRS & trs) {
		using std::distance;
		std::size_t vertexCount(0), edgeCount(0), maxVertexId(0);
		for(Finite_vertices_iterator it(trs.finite_vertices_begin()), end(trs.finite_vertices_end()); it != end; ++it) {
			if (!trs.is_special(it) && !trs.is_special(it)) {
				const VertexInfo & vi = it->info();
				assert(vi.valid());
				++vertexCount;
				maxVertexId = std::max<std::size_t>(maxVertexId, vi.id);
			}
		}
		for(Finite_edges_iterator it(trs.finite_edges_begin()), end(trs.finite_edges_end()) ; it != end; ++it) {
			auto & face = *(it->first);
			int index = it->second;
			auto vs = face.vertex(face.cw(index));
			auto vt = face.vertex(face.ccw(index));
			if (trs.is_special(vs) || trs.is_special(vt)) {
				continue;
			}
			++edgeCount;
		}
		
		out << vertexCount << '\n';
		out << edgeCount << '\n';
		
		ratss::OutputPoint op(3);
		std::vector<int> vertexId2OutId(maxVertexId+1);
		std::size_t counter = 0;
		for(Finite_vertices_iterator it(trs.finite_vertices_begin()), end(trs.finite_vertices_end()); it != end; ++it) {
			if (trs.is_special(it)) {
				continue;
			}
		
			const VertexInfo & vi = it->info();
			const Point & p = it->point();
			op.coords[0] = ratss::Conversion<FT>::toMpq( p.x() );
			op.coords[1] = ratss::Conversion<FT>::toMpq( p.y() );
			op.coords[2] = ratss::Conversion<FT>::toMpq( p.z() );
			
			op.print(out, outFormat);
			out << '\n';
			vertexId2OutId.at(vi.id) = counter;
			++counter;
		}
		
		for(Finite_edges_iterator it(trs.finite_edges_begin()), end(trs.finite_edges_end()) ; it != end; ++it) {
			auto & face = *(it->first);
			int index = it->second;
			auto vs = face.vertex(face.cw(index));
			auto vt = face.vertex(face.ccw(index));
			if (trs.is_special(vs) || trs.is_special(vt)) {
				continue;
			}
			
			const VertexInfo & vis = vs->info();
			const VertexInfo & vit = vt->info();
			out << vis.id << ' ' << vit.id << '\n';
		}
	}

	void writeSimplestGraphRendering(std::ostream & out, TRS & sdt) {
		out.precision(20);
		
		CGAL::Unique_hash_map<Vertex_handle, uint32_t> vertex2Id;
		uint32_t vertexId = 0;
		uint32_t faceCount = 0;
		
		Finite_vertices_iterator vIt(sdt.finite_vertices_begin()), vEnd(sdt.finite_vertices_end());
		for(; vIt != vEnd; ++vIt) {
			assert(!vertex2Id.is_defined(vIt));
			vertex2Id[vIt] = vertexId;
			++vertexId;
		}
		
		Finite_faces_iterator fIt(sdt.finite_faces_begin()), fEnd(sdt.finite_faces_end());
		for(; fIt != fEnd; ++fIt) {
			++faceCount;
		}
		
		std::cout << "Graph has " << vertexId << " vertices, " << faceCount << '/' << faceCount*3 << " faces/edges" << std::endl;
		
		out << vertexId << '\n';
		out << faceCount*3 << '\n';
		
		vIt = sdt.finite_vertices_begin();
		for(; vIt != vEnd; ++vIt) {
			assert(vertex2Id.is_defined(vIt));
			auto coords = sdt.toGeo(vIt->point());
			out << coords.lat << ' ' << coords.lon << '\n';
		}
		
		//and now the edges
		fIt = sdt.finite_faces_begin();
		for(; fIt != fEnd; ++fIt) {
			for(int j(0); j < 3; ++j) {
				auto vhs = fIt->vertex(j);
				auto vhe = fIt->vertex(sdt.ccw(j));
				assert(vertex2Id.is_defined(vhs));
				assert(vertex2Id.is_defined(vhe));
				
				out << vertex2Id[vhs] << ' ' << vertex2Id[vhe] << ' ' << 1 << ' ' << 1 << '\n';
			}
		}
		out << std::flush;
		
	}

	// #Knoten
	// #Kanten
	// Knotenliste: lat lon pro Zeile
	// Kantenliste: Start Ziel Breite Farbe pro Zeile
	// 
	// Farbe ist 1-5
	// start/ziel als index
	// lat/lon als float
	void writeSimplestGraphRenderingAndre(std::ostream & out, TRS & sdt) {
		out.precision(20);

// 		std::string node_rgb_value = "0 0 0 0";
// 		std::string triang_edge_rgb_value = "255 255 255 100";
// 		std::string constraint_edge_rgb_value = "255 255 255 255";
// 		std::string intersection_edge_rgb_value = "255 0 0 255";
// 		std::string triangle_rgb_value = "255 255 255 0";

		std::string node_rgb_value = "0 0 0 0";
		std::string triang_edge_rgb_value = "0 0 0 100";
		std::string constraint_edge_rgb_value = "0 0 0 255";
		std::string intersection_edge_rgb_value = "255 0 0 255";
		std::string triangle_rgb_value = "255 255 255 0";
		
		CGAL::Unique_hash_map<Vertex_handle, uint32_t> vertex2Id;
		uint32_t vertexId = 0;
		uint32_t edgeCount = 0;
		uint32_t faceCount = 0;
		
		Finite_vertices_iterator vIt(sdt.finite_vertices_begin()), vEnd(sdt.finite_vertices_end());
		for(; vIt != vEnd; ++vIt) {
			assert(!vertex2Id.is_defined(vIt));
			vertex2Id[vIt] = vertexId;
			++vertexId;
		}
		
		Finite_edges_iterator eIt(sdt.finite_edges_begin()), eEnd(sdt.finite_edges_end());
		for(; eIt != eEnd; ++eIt) {
			++edgeCount;
		}
		
		Finite_faces_iterator fIt(sdt.finite_faces_begin()), fEnd(sdt.finite_faces_end());
		for(; fIt != fEnd; ++fIt) {
			++faceCount;
		}
		
		std::cout << "Graph has " << vertexId << " vertices, " << faceCount << '/' << faceCount*3 << " faces/edges" << std::endl;
		
		out << vertexId << '\n';
		out << edgeCount << '\n';
		out << faceCount << '\n';
		
		vIt = sdt.finite_vertices_begin();
		for(; vIt != vEnd; ++vIt) {
			assert(vertex2Id.is_defined(vIt));
			auto coords = sdt.toGeo(vIt->point());
			out << coords.lat << ' ' << coords.lon << ' ' << node_rgb_value << '\n';
		}
		
		//and now the edges
		eIt = sdt.finite_edges_begin();
		for(; eIt != eEnd; ++eIt) {
			const Edge & e = *eIt;
			auto & face = *(eIt->first);
			int index = eIt->second;
			auto vs = face.vertex(face.cw(index));
			auto vt = face.vertex(face.ccw(index));
			assert(vertex2Id.is_defined(vs) && vertex2Id.is_defined(vt));
			out << vertex2Id[vs] << ' ' << vertex2Id[vt] << ' ';
			if (::is_constrained(sdt, e)) {
				if (!vs->info().valid() || !vt->info().valid()) {
					out << intersection_edge_rgb_value;
				}
				else {
					out << constraint_edge_rgb_value;
				}
			}
			else {
				out << triang_edge_rgb_value;
			}
			out << '\n';
		}
		
		//and now the faces
		fIt = sdt.finite_faces_begin();
		for(; fIt != fEnd; ++fIt) {
			auto v0 = fIt->vertex(0);
			auto v1 = fIt->vertex(1);
			auto v2 = fIt->vertex(2);
			
			assert(vertex2Id.is_defined(v0));
			assert(vertex2Id.is_defined(v1));
			assert(vertex2Id.is_defined(v2));
			
			out << vertex2Id[v0] << ' '
				<< vertex2Id[v1] << ' '
				<< vertex2Id[v2] << ' '
				<< triangle_rgb_value << '\n';
		}
		out << std::flush;
	}
	
};

class TriangulationCreator {
public:
	GraphOutputType got;
	ratss::OutputPoint::Format pointFormat;
	using Points = std::vector<std::pair<Point3, VertexInfo>>;
	using Edges = std::vector<std::pair<int, int>>;
public:
	TriangulationCreator() : got(GOT_INVALID), pointFormat(ratss::OutputPoint::FM_INVALID) {}
	virtual ~TriangulationCreator() {}
public:
	virtual void create(Points & points, Edges & edges, InputOutput & io, bool clear) = 0;
	virtual void write(InputOutput & io) = 0;
};

class TriangulationCreatorDelaunay: public TriangulationCreator {
public:
	using Tr =  dts2::Delaunay_triangulation_with_info_s2<VertexInfo, void>;
	using Point_3 = Tr::Point_3;
private:
	Tr m_tr;
public:
	TriangulationCreatorDelaunay(int significands) : m_tr(significands) {}
	
	virtual void create(Points & points, Edges & edges, InputOutput & io, bool clear) override {
		m_tr.insert(points.begin(), points.end());
		if (clear) {
			points = Points();
			edges = Edges();
		}
	}
	
	virtual void write(InputOutput & io) override {
		TriangulationWriter<Tr> writer(pointFormat, got);
		writer.write(io.output(), m_tr);
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
	
	virtual void create(Points & points, Edges & edges, InputOutput & io, bool clear) override  {
		std::size_t ps = points.size();
		m_tr.insert(points.begin(), points.end());
		if (clear) {
			points = Points();
		}
		
		std::vector<Vertex_handle> pId2Vertex(ps);
		for(Finite_vertices_iterator it(m_tr.finite_vertices_begin()), end(m_tr.finite_vertices_end()); it != end; ++it) {
			const VertexInfo & vi = it->info();
			if (vi.valid()) {
				Vertex_handle & vh = pId2Vertex.at(vi.id);
				vh =(Vertex_handle) it;
			}
		}
		Vertex_handle nullHandle;
		for(const std::pair<int, int> & e : edges) {
			if (e.first != e.second) {
				const auto & p1 = pId2Vertex.at(e.first);
				const auto & p2 = pId2Vertex.at(e.second);
				if (p1 != p2 && p1 != nullHandle && p2 != nullHandle) {
					m_tr.insert(p1, p2);
				}
			}
		}
		if (clear) {
			edges = Edges();
		}
	}
	
	virtual void write(InputOutput & io) override  {
		TriangulationWriter<Tr> writer(pointFormat, got);
		writer.write(io.output(), m_tr);
	}
};

using TriangulationCreatorNoIntersectionsConstrainedDelaunay = TriangulationCreatorConstrainedDelaunay<dts2::Constrained_Delaunay_triangulation_no_intersections_with_info_s2>;
using TriangulationCreatorInExactIntersectionsConstrainedDelaunay = TriangulationCreatorConstrainedDelaunay<dts2::Constrained_Delaunay_triangulation_with_inexact_intersections_with_info_s2>;

class TriangulationCreatorExactIntersectionsConstrainedDelaunay: public TriangulationCreator {
public:
	using Tr =  dts2::Constrained_Delaunay_triangulation_with_exact_intersections_with_info_s2<VertexInfo, void>;
	using Vertex_handle = typename Tr::Vertex_handle;
	using Finite_vertices_iterator = typename Tr::Finite_vertices_iterator;
	using Point_3 = Tr::Point_3;
private:
	Tr m_tr;
public:
	TriangulationCreatorExactIntersectionsConstrainedDelaunay(int significands) : m_tr(significands) {}
	
	virtual void create(Points & points, Edges & edges, InputOutput & io, bool clear) override {
		std::size_t ps = points.size();
		std::vector< std::pair<Point_3, VertexInfo> > myPoints;
		myPoints.reserve(ps);
		for(const std::pair<Point3, VertexInfo> & pi : points) {
			CORE::Expr x = ratss::Conversion<CORE::Expr>::moveFrom( ratss::Conversion<K::FT>::toMpq(pi.first.x()) );
			CORE::Expr y = ratss::Conversion<CORE::Expr>::moveFrom( ratss::Conversion<K::FT>::toMpq(pi.first.y()) );
			CORE::Expr z = ratss::Conversion<CORE::Expr>::moveFrom( ratss::Conversion<K::FT>::toMpq(pi.first.z()) );
			Point_3 myPoint(x, y, z);
			myPoints.emplace_back(std::move(myPoint), pi.second);
		}
		if (clear) {
			points = Points();
		}
		m_tr.insert(myPoints.begin(), myPoints.end());
		myPoints = std::vector< std::pair<Point_3, VertexInfo> >();
		
		std::vector<Vertex_handle> pId2Vertex(ps);
		for(Finite_vertices_iterator it(m_tr.finite_vertices_begin()), end(m_tr.finite_vertices_end()); it != end; ++it) {
			const VertexInfo & vi = it->info();
			if (vi.valid()) {
				Vertex_handle & vh = pId2Vertex.at(vi.id);
				vh =(Vertex_handle) it;
			}
		}
		Vertex_handle nullHandle;
		for(const std::pair<int, int> & e : edges) {
			if (e.first != e.second) {
				const auto & p1 = pId2Vertex.at(e.first);
				const auto & p2 = pId2Vertex.at(e.second);
				if (p1 != p2 && p1 != nullHandle && p2 != nullHandle) {
					m_tr.insert(p1, p2);
				}
			}
		}
		if (clear) {
			edges = Edges();
		}
	}
	
	virtual void write(InputOutput & io) override {
		TriangulationWriter<Tr> writer(pointFormat, got);
		writer.write(io.output(), m_tr);
	}
};


class TriangulationCreatorExactIntersectionsSphericalConstrainedDelaunay: public TriangulationCreator {
public:
	using Tr = dts2::Constrained_Delaunay_triangulation_with_inexact_intersections_with_info_s2<VertexInfo, void>;
// 	using Tr =  dts2::Constrained_Delaunay_triangulation_with_exact_intersections_spherical_with_info_s2<VertexInfo, void>;
	using Vertex_handle = typename Tr::Vertex_handle;
	using Finite_vertices_iterator = typename Tr::Finite_vertices_iterator;
	using Point_3 = typename Tr::Point_3;
	using FT = typename Tr::FT;
private:
	Tr m_tr;
public:
	TriangulationCreatorExactIntersectionsSphericalConstrainedDelaunay(int significands) : m_tr(significands) {}
	
	virtual void create(Points & points, Edges & edges, InputOutput & io, bool clear) override {
		std::size_t ps = points.size();
		m_tr.insert(points.begin(), points.end());
		if (clear) {
			points = Points();
		}
		
		std::vector<Vertex_handle> pId2Vertex(ps);
		for(Finite_vertices_iterator it(m_tr.finite_vertices_begin()), end(m_tr.finite_vertices_end()); it != end; ++it) {
			const VertexInfo & vi = it->info();
			if (vi.valid()) {
				Vertex_handle & vh = pId2Vertex.at(vi.id);
				vh =(Vertex_handle) it;
			}
		}
		Vertex_handle nullHandle;
		for(const std::pair<int, int> & e : edges) {
			if (e.first != e.second) {
				const auto & p1 = pId2Vertex.at(e.first);
				const auto & p2 = pId2Vertex.at(e.second);
				if (p1 != p2 && p1 != nullHandle && p2 != nullHandle) {
					m_tr.insert(p1, p2);
				}
			}
		}
		if (clear) {
			edges = Edges();
		}
	}
	
	virtual void write(InputOutput & io) override {
		TriangulationWriter<Tr> writer(pointFormat, got);
		writer.write(io.output(), m_tr);
	}
};

class Config: public ratss::BasicCmdLineOptions {
public:
	TriangulationType triangType;
	GraphOutputType got;
public:
	using ratss::BasicCmdLineOptions::parse;
public:
	virtual bool parse(const std::string & currentToken,int & i, int argc, char ** argv) override;
	void help(std::ostream & out) const;
	void print(std::ostream & out) const;
};

struct Data {
	Data();
	~Data();
	std::vector<std::pair<Point3, VertexInfo>> points;
	std::vector<std::pair<int, int>> edges;
	TriangulationCreator * tc;
	
	void init(const Config& cfg);
	void read(InputOutput& io, const Config& cfg);
	void create(InputOutput& io);
	void write(InputOutput& io, const Config& cfg);
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

	io.setInput(cfg.inFileName);
	io.setOutput(cfg.outFileName);
	
	if (cfg.verbose) {
		cfg.print(io.info());
		io.info() << std::endl;
	}
	
	
	io.info() << "Initializing triangulation..." << std::endl;
	data.init(cfg);
	
	io.info() << "Reading graph..." << std::endl;
	data.read(io, cfg);
	
	io.info() << "Creating triangulation..." << std::endl;
	data.create(io);
	
	io.info() << "Writing triangulation..." << std::endl;
	data.write(io, cfg);
	
	io.info() << "All operations succeeded. Exiting" << std::endl;
	return 0;
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
		else if (type == "cxs" || type == "constrained-intersection-exact-spherical") {
			triangType = TT_CONSTRAINED_EXACT_SPHERICAL;
		}
		else {
			throw ratss::BasicCmdLineOptions::ParseError("Unknown triangulation type: " + type);
		}
		++i;
	}
	else if (token == "-g" && i+1 < argc) {
		std::string type( argv[i+1] );
		if (type == "wx") {
			got = GOT_WITHOUT_SPECIAL;
		}
		else if (type == "simplest") {
			got = GOT_SIMPLEST_GRAPH_RENDERING;
		}
		else if (type == "simplest_andre") {
			got = GOT_SIMPLEST_GRAPH_RENDERING_ANDRE;
		}
		else {
			throw ratss::BasicCmdLineOptions::ParseError("Unknown graph output type: " + type);
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
		"\t-t type\ttype = [d,delaunay, c,constrained,cx,constrained-intersection,cxe,constrained-intesection-exact, cxs, constrained-intersection-exact-spherical]\n"
		"\t-g type\ttype = [wx, witout_special, simplest, simplest_andre]\n";
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
	case TT_CONSTRAINED_EXACT_SPHERICAL:
		out << "constrained exact intersection using spherical kernel";
		break;
	default:
		out << "invalid";
		break;
	}
	out << '\n';
	out << "Graph output type: ";
	switch (got) {
	case GOT_WITHOUT_SPECIAL:
		out << "without special";
		break;
	case GOT_SIMPLEST_GRAPH_RENDERING:
		out << "simplest graph rendering";
		break;
	case GOT_SIMPLEST_GRAPH_RENDERING_ANDRE:
		out << "simplest graph rendering andre";
		break;
	default:
		out << "invalid";
		break;
	};
	out << '\n';
	ratss::BasicCmdLineOptions::options_selection(out);
}

Data::Data() : tc(0) {}

Data::~Data() {
	delete tc;
	tc = 0;
}

void Data::init(const Config & cfg) {
	if (tc) {
		delete tc;
		tc = 0;
	}
	switch (cfg.triangType) {
	case TT_DELAUNAY:
		tc = new TriangulationCreatorDelaunay(cfg.significands);
		break;
	case TT_CONSTRAINED:
		tc = new TriangulationCreatorNoIntersectionsConstrainedDelaunay(cfg.significands);
		break;
	case TT_CONSTRAINED_INEXACT:
		tc = new TriangulationCreatorInExactIntersectionsConstrainedDelaunay(cfg.significands);
		break;
	case TT_CONSTRAINED_EXACT:
		tc = new TriangulationCreatorExactIntersectionsConstrainedDelaunay(cfg.significands);
		break;
	case TT_CONSTRAINED_EXACT_SPHERICAL:
		tc = new TriangulationCreatorExactIntersectionsSphericalConstrainedDelaunay(cfg.significands);
		break;
	default:
		throw std::runtime_error("Unkown triangulation type");
	};
}

void Data::read(InputOutput & io, const Config & cfg) {
	std::istream & is = io.input();

	std::size_t num_points, num_edges;
	is >> num_points >> num_edges;
	
	io.info() << "Need to fetch " << num_points << " points and " << num_edges << " edges" << std::endl;
	
	points.reserve(num_points);
	edges.reserve(num_edges);
	
	ratss::InputPoint ip;
	ratss::OutputPoint op;
	ratss::ProjectSN proj;
	
	if (cfg.progress) {
		io.info() << std::endl;
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
			io.info() << '\xd' << counter/1000 << "k" << std::flush;
		}
	}
	if (cfg.progress) {
		io.info() << std::endl;
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

void Data::create(InputOutput & io) {
	tc->create(points, edges, io, true);
	points.clear();
	edges.clear();
}

void Data::write(InputOutput & io, const Config & cfg) {
	tc->got = cfg.got;
	tc->pointFormat = cfg.outFormat;
	tc->write(io);
}