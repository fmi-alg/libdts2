
/**
  * Graph format is as follows
  * number_of_points
  * number_of_constraints
  * points
  * constraints
  *
  * Points may be encoded in any format accepted my ratss::FloatPoint
  * Constraints reference points: id1 id2
  * 
  */

#include <time.h>
#include <sys/time.h>
#include <fstream>
#include <string>
#include <libdts2/Constrained_delaunay_triangulation_s2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Spatial_sort_traits_adapter_3.h>
#include <libratss/mpreal.h>
#include <libratss/GeoCoord.h>
#include <libratss/SphericalCoord.h>
#include <libratss/ProjectSN.h>
#include <libratss/util/BasicCmdLineOptions.h>
#include <libratss/util/InputOutputPoints.h>
#include <libratss/util/InputOutput.h>
#include <libratss/Conversion.h>
#include <libratss/debug.h>
#ifndef __APPLE__
	#include <malloc.h>
#else
	void malloc_trim(int) {}
#endif

template<typename T_VERTEX_INFO, typename T_FACE_INFO>
using Delaunay_triangulation_with_info_s2_epeck = dts2::Delaunay_triangulation_with_info_s2<T_VERTEX_INFO, T_FACE_INFO, CGAL::Exact_predicates_exact_constructions_kernel>;

template<typename T_VERTEX_INFO, typename T_FACE_INFO>
using Delaunay_triangulation_with_info_s2_fsceik =
	dts2::Delaunay_triangulation_with_info_s2<T_VERTEX_INFO, T_FACE_INFO, CGAL::Filtered_simple_cartesian_extended_integer_kernel>;
	
template<typename T_VERTEX_INFO, typename T_FACE_INFO>
using Delaunay_triangulation_with_info_s2_spherical =
	dts2::Delaunay_triangulation_s2<
		dts2::Constrained_delaunay_triangulation_with_exact_intersections_spherical_traits_s2,
		typename dts2::internal::TriangulationDataStructureSelector<
			dts2::Constrained_delaunay_triangulation_with_exact_intersections_spherical_traits_s2,
			T_VERTEX_INFO,
			T_FACE_INFO
		>::type
	>;

class MemUsage {
public:
	MemUsage() {
		update();
	}
	void update() {
		const char* statm_path = "/proc/self/statm";

		FILE * f = fopen(statm_path,"r");
		if(!f){
			perror(statm_path);
			return;
		}
		if(7 != fscanf(f,"%ld %ld %ld %ld %ld %ld %ld", &size,&resident,&share,&text,&lib,&data,&dt)) {
			perror(statm_path);
		}
		fclose(f);
	}
	unsigned long size,resident,share,text,lib,data,dt;
	static unsigned long pagesize;
};

std::ostream & operator<<(std::ostream & out, const MemUsage & mem) {
	out << (mem.resident * mem.pagesize)/(1024*1024) << " MiB";
	return out;
}

class TimeMeasurer {
private:
	struct timeval m_begin, m_end;
public:
	TimeMeasurer() {
		::memset(&m_begin, 0, sizeof(struct timeval));
		::memset(&m_end, 0, sizeof(struct timeval));
	}
	
	~TimeMeasurer() {}
	
	inline void begin() {
		gettimeofday(&m_begin, NULL);
	}
	
	inline void end() {
		gettimeofday(&m_end, NULL);
	}
	
	inline long beginTime() const {
		return m_begin.tv_sec;
	}
	
	/** @return returns the elapsed time in useconds  */
	inline long elapsedTime() const {
		long mtime, seconds, useconds;
		seconds  = m_end.tv_sec  - m_begin.tv_sec;
		useconds = m_end.tv_usec - m_begin.tv_usec;
		mtime = (long)((double)((seconds) * 1000*1000 + useconds) + 0.5);
		return mtime;
	}
	
	inline long elapsedUseconds() const {
		return elapsedTime();
	}

	inline long elapsedMilliSeconds() const {
		return elapsedTime()/1000;
	}
	
	inline long elapsedSeconds() const {
		return elapsedTime()/1000000;
	}

	inline long elapsedMinutes() const {
		return elapsedSeconds()/60;
	}
	inline long elapsedHours() const {
		return elapsedSeconds()/3600;
	}
};

std::ostream & operator<<(std::ostream & out, const TimeMeasurer & tm) {
	long tms[4];
	tms[0] = tm.elapsedMilliSeconds()%1000;
	tms[1] = tm.elapsedSeconds()%60; //seconds
	tms[2] = tm.elapsedMinutes()%60; //minutes
	tms[3] = tm.elapsedHours(); //hours
	std::string e[4] = {"ms", "s", "M", "h"};
	
	bool hasPrev = false;
	for(int i = 3; i >= 0; --i) {
		if (tms[i]) {
			if (hasPrev)
				out << " ";
			hasPrev = true;
			out << tms[i] << e[i];
		}
	}
	if (!hasPrev) {
		out << "0s";
	}
	return out;
}

unsigned long MemUsage::pagesize = ::sysconf(_SC_PAGE_SIZE);

class Point3 {
public:
	using Epick = CGAL::Exact_predicates_inexact_constructions_kernel;
	using EpickPoint = Epick::Point_3;
public:
	using Epeck = CGAL::Exact_predicates_exact_constructions_kernel;
	using EpeckPoint = Epeck::Point_3;
public:
	using Sceik = CGAL::Simple_cartesian_extended_integer_kernel;
	using Fsceik = CGAL::Filtered_simple_cartesian_extended_integer_kernel;
	using Flceik = CGAL::Filtered_lazy_cartesian_extended_integer_kernel;
	using Epeceik = CGAL::Exact_predicates_exact_constructions_extended_integer_kernel;
	
	using SceikPoint = Sceik::Point_3;
	using FsceikPoint = Fsceik::Point_3;
	using FlceikPoint = Flceik::Point_3;
	using EpeceikPoint = Epeceik::Point_3;
public:
	using Sce1024ik = CGAL::Simple_cartesian_extended_1024_integer_kernel;
	using Fsce1024ik = CGAL::Filtered_simple_cartesian_extended_1024_integer_kernel;
	using Flce1024ik = CGAL::Filtered_lazy_cartesian_extended_1024_integer_kernel;
	
	using Fsce1024ikPoint = Fsce1024ik::Point_3;
	using Flce1024ikPoint = Flce1024ik::Point_3;
public:
	using Epecksqrt = CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt;
	using EpecksqrtPoint = Epecksqrt::Point_3;
public:
	using Epecsk = CGAL::Exact_spherical_kernel_3;
	using EpecskPoint = Epecsk::Circular_arc_point_3;
	using EpecslkPoint = Epecsk::Linear_kernel::Point_3;
public:
// 	using FT = mpq_class;
	using K = Sceik;
	using FT = K::FT;
	using R = K::Point_3::R;
public:
	Point3(const FT & x, const FT & y, const FT & z) :
	m_x(x), m_y(y), m_z(z)
	{}
	template<typename T_POINT>
	Point3(const T_POINT & p) :
	m_x( ratss::convert<FT>( p.x() )),
	m_y( ratss::convert<FT>( p.y() )),
	m_z( ratss::convert<FT>( p.z() ))
	{}
	~Point3() {}
public:
	const FT & x() const { return m_x; }
	const FT & y() const { return m_y; }
	const FT & z() const { return m_z; }
public:
	template<typename T_KERNEL>
	typename T_KERNEL::Point_3 convert_to_p() const {
		return typename T_KERNEL::Point_3(
			ratss::convert<typename T_KERNEL::FT>( x() ),
			ratss::convert<typename T_KERNEL::FT>( y() ),
			ratss::convert<typename T_KERNEL::FT>( z() )
		);
	}
public:
	explicit operator EpickPoint() const { return convert_to_p<Epick>(); }
	explicit operator EpeckPoint() const { return convert_to_p<Epeck>(); }
	explicit operator SceikPoint() const { return convert_to_p<Sceik>(); }
	explicit operator FsceikPoint() const { return convert_to_p<Fsceik>(); }
	explicit operator FlceikPoint() const { return convert_to_p<Flceik>(); }
	explicit operator EpeceikPoint() const { return convert_to_p<Epeceik>(); }
	explicit operator Fsce1024ikPoint() const { return convert_to_p<Fsce1024ik>(); }
	explicit operator Flce1024ikPoint() const { return convert_to_p<Flce1024ik>(); }
	explicit operator EpecksqrtPoint() const { return convert_to_p<Epecksqrt>(); }
	explicit operator EpecslkPoint() const { return convert_to_p<Epecsk::Linear_kernel>(); } 
	explicit operator EpecskPoint() const { return EpecskPoint( static_cast<EpecslkPoint>(*this) );}
	template<typename T_KERNEL>
	explicit operator dts2::Point_sp<T_KERNEL> () const {
		return dts2::Point_sp<T_KERNEL>( convert_to_p<T_KERNEL>() );
	}
private:
	FT m_x;
	FT m_y;
	FT m_z;
	
};

class BinaryIo {
public:
	BinaryIo() : m_io(0) {}
	void setIo(ratss::InputOutput * io) {
		m_io = io;
	}
public:
	template<typename T>
	T get() { return get(static_cast<T*>(0)); };
	template<typename T>
	void write(const T & v);
private:
	template<typename T>
	T get(T *); //we need this, since we cannot simply specialize this based on T, since it is only the return type
	CGAL::ExtendedInt64q<CGAL::Gmpq>::PQ get(CGAL::ExtendedInt64q<CGAL::Gmpq>::PQ*);
	CGAL::ExtendedInt64q<CGAL::Gmpq> get(CGAL::ExtendedInt64q<CGAL::Gmpq>*);
	Point3 get(Point3*);
private:
	ratss::InputOutput * m_io;
};

template<typename T>
T
BinaryIo::get(T*) {
	assert(std::is_pod<T>::value);
	char tmp[sizeof(T)];
	m_io->input().read(tmp, sizeof(T));
	T result;
	::memmove(&result, tmp, sizeof(T));
	return result;
}

template<typename T>
void
BinaryIo::write(T const & v) {
	assert(std::is_pod<T>::value);
	char tmp[sizeof(T)];
	::memmove(tmp, &v, sizeof(T));
	m_io->output().write(tmp, sizeof(T));
}

CGAL::ExtendedInt64q<CGAL::Gmpq>::PQ
BinaryIo::get(CGAL::ExtendedInt64q<CGAL::Gmpq>::PQ* = 0) {
	int64_t num = this->get<int64_t>();
	int64_t den = this->get<int64_t>();
	return CGAL::ExtendedInt64q<CGAL::Gmpq>::PQ(num, den);
}

CGAL::ExtendedInt64q<CGAL::Gmpq>
BinaryIo::get(CGAL::ExtendedInt64q<CGAL::Gmpq>* = 0) {
	return CGAL::ExtendedInt64q<CGAL::Gmpq>( this->get<CGAL::ExtendedInt64q<CGAL::Gmpq>::PQ>() );
}

Point3
BinaryIo::get(Point3* = 0) {
	using FT = CGAL::ExtendedInt64q<CGAL::Gmpq>;
	auto x = this->get<FT>();
	auto y = this->get<FT>();
	auto z = this->get<FT>();
	return Point3(
		ratss::convert<Point3::FT>( x ),
		ratss::convert<Point3::FT>( y ),
		ratss::convert<Point3::FT>( z )
	);
}

template<>
void
BinaryIo::write(const CGAL::ExtendedInt64q<CGAL::Gmpq>::PQ & v) {
	this->write<int64_t>(v.num);
	this->write<int64_t>(v.den);
}

template<>
void
BinaryIo::write(const CGAL::ExtendedInt64q<CGAL::Gmpq> & v) {
	if (v.isExtended()) {
		throw std::domain_error("Unable to serialize extended CGAL::ExtendedInt64q<CGAL::Gmpq>");
	}
	this->write(v.getPq());
}

template<>
void
BinaryIo::write(const Point3 & v) {
	using FT = CGAL::ExtendedInt64q<CGAL::Gmpq>;
	this->write( ratss::convert<FT>( v.x() ) );
	this->write( ratss::convert<FT>( v.y() ) );
	this->write( ratss::convert<FT>( v.z() ) );
}

typedef enum {TT_DELAUNAY, TT_DELAUNAY_64, TT_DELAUNAY_SPHERICAL, TT_CONVEX_HULL_INEXACT, TT_CONVEX_HULL, TT_CONVEX_HULL_64, TT_CONSTRAINED, TT_CONSTRAINED_INEXACT, TT_CONSTRAINED_INEXACT_64, TT_CONSTRAINED_INEXACT_SP, TT_CONSTRAINED_INEXACT_SPK, TT_CONSTRAINED_INEXACT_SPK64, TT_CONSTRAINED_EXACT, TT_CONSTRAINED_EXACT_SPHERICAL} TriangulationType;
typedef enum {GOT_INVALID, GOT_NONE, GOT_WITHOUT_SPECIAL, GOT_WITHOUT_SPECIAL_HASH_MAP, GOT_SIMPLEST_GRAPH_RENDERING, GOT_SIMPLEST_GRAPH_RENDERING_ANDRE} GraphOutputType;
typedef enum {GIT_INVALID, GIT_NODES_EDGES, GIT_EDGES, GIT_NODES_EDGES_BY_POINTS_BINARY} GraphInputType;
typedef enum {TIO_INVALID, TIO_NODES_EDGES, TIO_EDGES, TIO_NODES_EDGES_BY_POINTS_BINARY} TriangulationInputOrder;

struct VertexInfo {
	VertexInfo() : id(def_instance_counter) { --def_instance_counter;}
	VertexInfo(int id) : id(id) {}
	bool valid() const { return id >= 0; }
	int id;
	static int def_instance_counter;
};

int VertexInfo::def_instance_counter = -1;

template<typename T_LINEAR_KERNEL>
bool
is_constrained(
	const dts2::Delaunay_triangulation_with_info_s2<VertexInfo, void, T_LINEAR_KERNEL> & /*trs*/,
	const typename dts2::Delaunay_triangulation_with_info_s2<VertexInfo, void, T_LINEAR_KERNEL>::Edge & /*e*/)
{
	return false;
}

bool
is_constrained(
	const Delaunay_triangulation_with_info_s2_spherical<VertexInfo, void> & /*trs*/,
	const Delaunay_triangulation_with_info_s2_spherical<VertexInfo, void>::Edge & /*e*/)
{
	return false;
}

template<typename T_TRS>
bool is_constrained(const T_TRS & trs, const typename T_TRS::Edge & e) {
	return trs.is_constrained(e);
}

extern "C" void debug_print_gmpq(const CGAL::Gmpq & pq) {
	std::cerr << pq << std::endl;
}

extern "C" void debug_print_point3(const Point3 &p ) {
	ratss::ProjectS2 proj;
	ratss::GeoCoord gc;
	proj.toGeo(p.x(), p.y(), p.z(), gc.lat, gc.lon, 53);
	std::cerr << p.x() << " " << p.y() << " " << p.z() << std::endl;
	std::streamsize prec = std::cerr.precision();
	std::cerr.precision(std::numeric_limits<double>::digits10+1);
	std::cerr << gc << std::endl;
	std::cerr.precision(prec);
}

extern "C" void debug_print_vertex_handle_icdt(const dts2::Constrained_Delaunay_triangulation_with_inexact_intersections_with_info_s2<VertexInfo, void>::Vertex_handle & vh) {
	if (vh == dts2::Constrained_Delaunay_triangulation_with_inexact_intersections_with_info_s2<VertexInfo, void>::Vertex_handle()) {
		std::cerr << "VertexHandle is null" << std::endl;
	}
	else {
		std::cerr << "VertexInfo=" << vh->info().id << std::endl;
		debug_print_point3(vh->point());
	}
}

dts2::Constrained_Delaunay_triangulation_with_inexact_intersections_with_info_s2<VertexInfo, void> * debug_triangulation_pointer = 0;

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
	
	ratss::RationalPoint::Format outFormat;
	GraphOutputType got;
	
	TriangulationWriter(ratss::RationalPoint::Format outFormat, GraphOutputType got) : outFormat(outFormat), got(got) {}
	
	void write(InputOutput & out, TRS & trs) {
		switch (got) {
		case GOT_NONE:
			break;
		case GOT_WITHOUT_SPECIAL:
			writeWithoutSpecial(out, trs);
			break;
		case GOT_WITHOUT_SPECIAL_HASH_MAP:
			writeWithoutSpecialHashMap(out, trs);
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
	
	void assignVertexIds(TRS & trs, std::vector<Vertex_handle> & id2Vh) {
		int count = 0;
		for(Finite_vertices_iterator it(trs.finite_vertices_begin()), end(trs.finite_vertices_end()); it != end; ++it, ++count) {
			it->info().id = count;
		}
		id2Vh.reserve(count);
		for(Finite_vertices_iterator it(trs.finite_vertices_begin()), end(trs.finite_vertices_end()); it != end; ++it, ++count) {
			id2Vh.at(it->info().id) = it;
		}
	}
	
	void writeWithoutSpecial(InputOutput & io, TRS & trs) {
		std::ostream & out = io.output();
		using std::distance;
// 		assignVertexIds(trs);
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
		
		ratss::RationalPoint op(3);
		std::vector<int> vertexId2OutId(maxVertexId+1);
		std::size_t counter = 0;
		for(Finite_vertices_iterator it(trs.finite_vertices_begin()), end(trs.finite_vertices_end()); it != end; ++it) {
			if (trs.is_special(it)) {
				continue;
			}
		
			const VertexInfo & vi = it->info();
			const Point & p = it->point();
			op.coords[0] = ratss::convert<mpq_class>( p.x() );
			op.coords[1] = ratss::convert<mpq_class>( p.y() );
			op.coords[2] = ratss::convert<mpq_class>( p.z() );
			
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

	void writeWithoutSpecialHashMap(InputOutput & io, TRS & sdt) {
		std::ostream & out = io.output();
		
		out.precision(std::numeric_limits<double>::digits10+1);
		
		CGAL::Unique_hash_map<Vertex_handle, uint32_t> vertex2Id;
		uint32_t vertexId = 0;
		uint32_t edgeCount = 0;
		
		Finite_vertices_iterator vIt(sdt.finite_vertices_begin()), vEnd(sdt.finite_vertices_end());
		for(; vIt != vEnd; ++vIt) {
			assert(!vertex2Id.is_defined(vIt));
			if (!sdt.is_auxiliary(vIt)) {
				vertex2Id[vIt] = vertexId;
				++vertexId;
			}
		}

		for(Finite_edges_iterator eIt(sdt.finite_edges_begin()), eEnd(sdt.finite_edges_end()); eIt != eEnd; ++eIt) {
			auto & face = *(eIt->first);
			int index = eIt->second;
			auto vs = face.vertex(face.cw(index));
			auto vt = face.vertex(face.ccw(index));
			if (vertex2Id.is_defined(vs) && vertex2Id.is_defined(vt)) { //src or target is aux/inf 
				++edgeCount;
			}
		}
		
		io.info() << "Graph has " << vertexId << " vertices and " << edgeCount  << " edges" << std::endl;
		
		out << vertexId << '\n';
		out << edgeCount << '\n';
		
		vIt = sdt.finite_vertices_begin();
		for(; vIt != vEnd; ++vIt) {
			if (vertex2Id.is_defined(vIt)) {
				auto coords = sdt.toGeo(vIt->point());
				out << coords.lat << ' ' << coords.lon << '\n';
			}
		}
		
		//and now the edges
		for(Finite_edges_iterator eIt(sdt.finite_edges_begin()), eEnd(sdt.finite_edges_end()); eIt != eEnd; ++eIt) {
			const Edge & e = *eIt;
			auto & face = *(eIt->first);
			int index = eIt->second;
			auto vs = face.vertex(face.cw(index));
			auto vt = face.vertex(face.ccw(index));
			if (!vertex2Id.is_defined(vs) || !vertex2Id.is_defined(vt)) { //src or target is aux/inf 
				continue;
			}
			out << vertex2Id[vs] << ' ' << vertex2Id[vt] << ' ';
			if (::is_constrained(sdt, e)) {
				out << "c";
			}
			else {
				out << "n";
			}
			out << '\n';
		}
		out << std::flush;
		
	}

	void writeSimplestGraphRendering(InputOutput & io, TRS & sdt) {
		std::ostream & out = io.output();
		out.precision(std::numeric_limits<double>::digits10+1);
		
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
		
		io.info() << "Graph has " << vertexId << " vertices, " << faceCount << '/' << faceCount*3 << " faces/edges" << std::endl;
		
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
	void writeSimplestGraphRenderingAndre(InputOutput & io, TRS & sdt) {
		std::ostream & out = io.output();
		out.precision(20);

// 		std::string node_rgb_value = "0 0 0 0";
// 		std::string triang_edge_rgb_value = "255 255 255 100";
// 		std::string constraint_edge_rgb_value = "255 255 255 255";
// 		std::string intersection_edge_rgb_value = "255 0 0 255";
// 		std::string triangle_rgb_value = "255 255 255 0";

		std::string node_rgb_value = "0 0 0 0";
		std::string triang_edge_rgb_value = "120 120 0 255";
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
		
		io.info() << "Graph has " << vertexId << " vertices, " << faceCount << '/' << faceCount*3 << " faces/edges" << std::endl;
		
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

extern "C" void debug_dump_triangulation() {
	if (!debug_triangulation_pointer) {
		std::cerr << "Pointer is null" << std::endl;
		return;
	}
	using CDT = dts2::Constrained_Delaunay_triangulation_with_inexact_intersections_with_info_s2<VertexInfo, void>;
	CDT & cdt = *debug_triangulation_pointer;
	TriangulationWriter<CDT> writer(ratss::PointBase::FM_GEO, GOT_WITHOUT_SPECIAL_HASH_MAP);
	std::cerr << "DEBUG_GRAPH_BEGIN" << std::endl;
	InputOutput io(std::cin, std::cerr);
	writer.write(io, cdt);
	std::cerr << "DEBUG_GRAPH_END" << std::endl;
}

class TriangulationCreator {
public:
	GraphOutputType got;
	ratss::RationalPoint::Format pointFormat;
	using Points = std::vector<std::pair<Point3, VertexInfo>>;
	using Edges = std::vector<std::pair<int, int>>;
public:
	TriangulationCreator() : got(GOT_INVALID), pointFormat(ratss::RationalPoint::FM_INVALID) {}
	virtual ~TriangulationCreator() {}
public:
	///@param points need to be exactly on the sphere
	virtual void create(Points & points, Edges & edges, InputOutput & io, bool clear) = 0;
	virtual void add(const Point3 & p) = 0;
	virtual void add(const Point3 & p1, const Point3 & p2) = 0;
	virtual void info(InputOutput &) = 0;
	virtual void write(InputOutput & io) = 0;
};

template<typename T_KERNEL>
class ConvexHullTriangulationCreator: public TriangulationCreator {
public:
	using K = T_KERNEL;
	using Point_3 = typename K::Point_3;
	using Tr = CGAL::Polyhedron_3<K>;
public:
	ConvexHullTriangulationCreator() {}
	virtual ~ConvexHullTriangulationCreator() {}
public:
	///@param points need to be exactly on the sphere
	virtual void create(Points & points, Edges & edges, InputOutput & io, bool clear) {
		auto tf = [](const Points::value_type & pi) {
			return (Point_3) pi.first;
		};
		using MyIterator = boost::transform_iterator<decltype(tf), Points::const_iterator>;
		
		if (clear) {
			edges = Edges();
		}
		io.info() << "Creating triangulation using convex hull algorithm..." << std::flush;
		CGAL::convex_hull_3(MyIterator(points.begin(), tf), MyIterator(points.end(), tf), m_tr);
		io.info() << "done" << std::endl;
		if (clear) {
			points = Points();
		}
	}
	virtual void add(const Point3 & /*p*/) {
		throw std::runtime_error("Convex hull based triangulation does not allow a dynamic creation");
	}
	virtual void add(const Point3 & /*p1*/, const Point3 & /*p2*/) {
		throw std::runtime_error("Convex hull based triangulation does not allow a dynamic creation");
	}
	virtual void write(InputOutput & /*io*/) {
		if (got != GOT_NONE) {
			throw std::runtime_error("Convex hull based triangulation does not support writing of triangulation");
		}
	}
	
	virtual void info(InputOutput & io) override {
		io.info() << "Convex hull has " << m_tr.size_of_vertices() << " vertices and " << m_tr.size_of_facets() << " facets" << std::endl;
	}
	
private:
	Tr m_tr;
};

using EpickConvexHullTriangulationCreator = ConvexHullTriangulationCreator<CGAL::Exact_predicates_inexact_constructions_kernel>;
using EpeckConvexHullTriangulationCreator = ConvexHullTriangulationCreator<CGAL::Exact_predicates_exact_constructions_kernel>;
using Ei64ConvexHullTriangulationCreator = ConvexHullTriangulationCreator<CGAL::Exact_predicates_exact_constructions_extended_integer_kernel>;

template< template<typename, typename> class T_TRS>
class TriangulationCreatorDelaunay: public TriangulationCreator {
public:
	using Tr =  T_TRS<VertexInfo, void>;
	using FT = typename Tr::FT;
	using Point_3 = typename Tr::Point_3;
	using Vertex_handle = typename Tr::Vertex_handle;
	using Finite_vertices_iterator = typename Tr::Finite_vertices_iterator;
	using Face_handle = typename Tr::Face_handle;
public:
	TriangulationCreatorDelaunay(int significands) : m_tr(significands) {}
	
	virtual void create(Points & points, Edges & edges, InputOutput & io, bool clear) override {
		if (clear) {
			edges = Edges();
		}
	
		auto tf = [](const Points::value_type & pi) -> std::pair<Point_3, VertexInfo> {
			return std::pair<Point_3, VertexInfo>(static_cast<Point_3>(pi.first), pi.second);
		};
		using MyIterator = boost::transform_iterator<decltype(tf), Points::const_iterator>;
		io.info() << "Inserting points..." << std::flush;
		m_tr.insert(MyIterator(points.begin(), tf), MyIterator(points.end(), tf), false);
		io.info() << "done" << std::endl;
		
		if (clear) {
			points = Points();
		}
	}

	virtual void add(const Point3 & p) override {
		m_lastAddVh = m_tr.insert( static_cast<Point_3>(p), (m_lastAddVh != Vertex_handle() ? m_lastAddVh->face() : Face_handle()), false);
	}
	
	virtual void add(const Point3 & p1, const Point3 & p2) override {
		add(p1);
		add(p2);
	}
	
	virtual void info(InputOutput & io) override {
		io.info() << "Triangulation has " << m_tr.number_of_vertices() << " vertices and " << m_tr.number_of_faces() << " faces" << std::endl;
	}

	virtual void write(InputOutput & io) override {
		TriangulationWriter<Tr> writer(pointFormat, got);
		writer.write(io, m_tr);
	}
private:
	Tr m_tr;
	Vertex_handle m_lastAddVh;
};

using TriangulationCreatorDelaunayEpeck =
	TriangulationCreatorDelaunay<Delaunay_triangulation_with_info_s2_epeck>;
	
using TriangulationCreatorDelaunayFsceik =
	TriangulationCreatorDelaunay<Delaunay_triangulation_with_info_s2_fsceik>;
	
using TriangulationCreatorDelaunaySpherical =
	TriangulationCreatorDelaunay<Delaunay_triangulation_with_info_s2_spherical>;

template< template<typename, typename> class T_TRS>
class TriangulationCreatorConstrainedDelaunay: public TriangulationCreator {
public:
	using Tr =  T_TRS<VertexInfo, void>;
	using FT = typename Tr::FT;
	using Point = typename Tr::Point;
	using Vertex_handle = typename Tr::Vertex_handle;
	using Finite_vertices_iterator = typename Tr::Finite_vertices_iterator;
	using Face_handle = typename Tr::Face_handle;
public:
	TriangulationCreatorConstrainedDelaunay(int significands) :
	m_tr(significands)
	{}
	TriangulationCreatorConstrainedDelaunay(int significands, int intersectionSignificands);
	
	~TriangulationCreatorConstrainedDelaunay() {}
	
	virtual void create(Points & points, Edges & edges, InputOutput & io, bool clear) override  {
		std::size_t ps = points.size();
		io.info() << "Inserting points..." << std::flush;
		TimeMeasurer tm;
		tm.begin();
		this->insert(points.begin(), points.end());
		tm.end();
		io.info() << tm << std::endl;
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
		io.info() << "Filtering constraints..." << std::flush;
		tm.begin();
		auto relevantEdges = compute_relevant_edges(points, pId2Vertex, edges, io);
		tm.end();
		io.info() << tm << std::endl;
		
		io.info() << "Inserting constraints..." << std::flush;
		tm.begin();
		insert_constraints(pId2Vertex, edges, relevantEdges, io);
		tm.end();
		io.info() << tm << std::endl;
		if (clear) {
			edges = Edges();
		}
	}
	std::vector<bool> compute_relevant_edges(Points & points, std::vector<Vertex_handle> & pId2Vertex, Edges & edges, InputOutput & /*io*/) {
		std::vector<bool> result(edges.size(), false);
		Vertex_handle nullHandle;
		Point3::FT maxLen(0.5);
		for(std::size_t i(0), s(edges.size()); i < s; ++i) {
			const std::pair<int, int> & e = edges[s-i-1];
			if (e.first != e.second) {
				const auto & v1 = pId2Vertex.at(e.first);
				const auto & v2 = pId2Vertex.at(e.second);
				if (v1 != v2 && v1 != nullHandle && v2 != nullHandle) {
					if (csd2(points.at(e.first).first, points.at(e.second).first) < maxLen) {
						result[i] = true;
					}
					else {
						std::cerr << "Removed long edge" << std::endl;
					}
				}
			}
		}
		return result;
	}
	
	auto csd2 (const Point3 & a, const Point3 & b) {
		auto x = a.x() - b.x();
		auto y = a.y() - b.y();
		auto z = a.z() - b.z();
		return x*x + y*y + z*z;
	}
	
	void insert_constraints(std::vector<Vertex_handle> & pId2Vertex, Edges & edges, std::vector<bool> const & relevantEdges, InputOutput & /*io*/) {
		Vertex_handle nullHandle;
		FT maxLen(0.5);
		for(std::size_t i(0), s(edges.size()); i < s; ++i) {
			if (relevantEdges[i]) {
				const std::pair<int, int> & e = edges[s-i-1];
				const auto & v1 = pId2Vertex.at(e.first);
				const auto & v2 = pId2Vertex.at(e.second);
				m_tr.insert(v1, v2);
			}
		}
	}
	
	virtual void add(const Point3 & p) override {
		m_lastAddVh = insert(p, (m_lastAddVh != Vertex_handle() ? m_lastAddVh->face() : Face_handle()));
	}
	
	virtual void add(const Point3 & p1, const Point3 & p2) override {
		add(p1);
		Vertex_handle vh1 = m_lastAddVh;
		add(p2);
		Vertex_handle vh2 = m_lastAddVh;
		insert_constraint(vh1, vh2);
		m_lastAddVh = vh1;
	}
	
	virtual void info(InputOutput & io) override {
		io.info() << "Triangulation has " << m_tr.number_of_vertices() << " vertices and " << m_tr.number_of_faces() << " faces" << std::endl;
	}
	
	virtual void write(InputOutput & io) override  {
		TriangulationWriter<Tr> writer(pointFormat, got);
		writer.write(io, m_tr);
	}
private:
	//points are already snapped
	void insert(Points::const_iterator begin, Points::const_iterator end) {
		auto tf = [](const Points::value_type & pi) {
			return std::pair<Point, VertexInfo>((Point) pi.first, pi.second);
		};
		using MyIterator = boost::transform_iterator<decltype(tf), Points::const_iterator>;
		m_tr.insert(MyIterator(begin, tf), MyIterator(end, tf), false);
	}
	Vertex_handle insert(const Point3 & p, const Face_handle & fh) {
		return m_tr.insert((Point) p, fh, false);
	}
	void insert_constraint(const Point3 & p1, const Point3 & p2) {
		m_tr.insert((Point) p1, (Point) p2);
	}
	void insert_constraint(const Vertex_handle & vh1, const Vertex_handle & vh2) {
		if (vh1 != vh2) {
			m_tr.insert(vh1, vh2);
		}
	}
private:
	Tr m_tr;
	Vertex_handle m_lastAddVh;
};

using TriangulationCreatorNoIntersectionsConstrainedDelaunay =
	TriangulationCreatorConstrainedDelaunay<dts2::Constrained_Delaunay_triangulation_no_intersections_with_info_s2>;
using TriangulationCreatorInExactIntersectionsConstrainedDelaunay =
	TriangulationCreatorConstrainedDelaunay<dts2::Constrained_Delaunay_triangulation_with_inexact_intersections_with_info_s2>;
using TriangulationCreatorInExactIntersectionsConstrainedDelaunay64 =
	TriangulationCreatorConstrainedDelaunay<dts2::Constrained_Delaunay_triangulation_with_inexact_intersections_with_info_s2_64>;
using TriangulationCreatorInExactIntersectionsConstrainedDelaunaySp =
	TriangulationCreatorConstrainedDelaunay<dts2::Constrained_Delaunay_triangulation_with_inexact_intersections_with_info_s2_sp>;
using TriangulationCreatorInExactIntersectionsConstrainedDelaunaySpk64 =
	TriangulationCreatorConstrainedDelaunay<dts2::Constrained_Delaunay_triangulation_with_inexact_intersections_with_info_s2_spk_64>;
using TriangulationCreatorInExactIntersectionsConstrainedDelaunaySpk =
	TriangulationCreatorConstrainedDelaunay<dts2::Constrained_Delaunay_triangulation_with_inexact_intersections_with_info_s2_spk>;
using TriangulationCreatorExactIntersectionsConstrainedDelaunay = 
	TriangulationCreatorConstrainedDelaunay<dts2::Constrained_Delaunay_triangulation_with_exact_intersections_with_info_s2>;

using TriangulationCreatorExactIntersectionsSphericalConstrainedDelaunay =
	TriangulationCreatorConstrainedDelaunay<dts2::Constrained_Delaunay_triangulation_with_exact_intersections_spherical_with_info_s2>;

// using TriangulationCreatorExactIntersectionsSphericalConstrainedDelaunay =
// 	TriangulationCreatorExactIntersectionsConstrainedDelaunay;
	
template<>
TriangulationCreatorNoIntersectionsConstrainedDelaunay::TriangulationCreatorConstrainedDelaunay(int significands, int /*intersectionSignificands*/) :
TriangulationCreatorConstrainedDelaunay(significands)
{}

template<>
TriangulationCreatorInExactIntersectionsConstrainedDelaunay::TriangulationCreatorConstrainedDelaunay(int significands, int intersectionSignificands) :
m_tr(
	TriangulationCreatorInExactIntersectionsConstrainedDelaunay::Tr::Geom_traits(
		LIB_RATSS_NAMESPACE::Conversion<TriangulationCreatorInExactIntersectionsConstrainedDelaunay::Tr::Geom_traits::FT>::moveFrom(
			LIB_RATSS_NAMESPACE::Conversion<uint64_t>::toMpq(std::numeric_limits<uint64_t>::max()-1)/
			LIB_RATSS_NAMESPACE::Conversion<uint64_t>::toMpq(std::numeric_limits<uint64_t>::max())
		),
		significands,
		intersectionSignificands
	)
)
{
	debug_triangulation_pointer = &m_tr;
}

template<>
TriangulationCreatorInExactIntersectionsConstrainedDelaunay64::TriangulationCreatorConstrainedDelaunay(int significands, int intersectionSignificands) :
m_tr(
	TriangulationCreatorInExactIntersectionsConstrainedDelaunay64::Tr::Geom_traits(
		LIB_RATSS_NAMESPACE::Conversion<TriangulationCreatorInExactIntersectionsConstrainedDelaunay64::Tr::Geom_traits::FT>::moveFrom(
			LIB_RATSS_NAMESPACE::Conversion<uint64_t>::toMpq(std::numeric_limits<uint64_t>::max()-1)/
			LIB_RATSS_NAMESPACE::Conversion<uint64_t>::toMpq(std::numeric_limits<uint64_t>::max())
		),
		significands,
		intersectionSignificands
	)
)
{}

template<>
TriangulationCreatorInExactIntersectionsConstrainedDelaunay64::~TriangulationCreatorConstrainedDelaunay() {
#ifdef WITH_EI64_STATS
	std::cerr << "ExtendedInt64q::number_of_allocations=" <<
		CGAL::ExtendedInt64q<CGAL::Gmpq>::number_of_allocations << std::endl;
	std::cerr << "ExtendedInt64q::number_of_extended_allocations=" <<
		CGAL::ExtendedInt64q<CGAL::Gmpq>::number_of_extended_allocations << std::endl;
	std::cerr << "ExtendedInt64z::number_of_extended_allocations=" <<
		CGAL::ExtendedInt64z::number_of_extended_allocations << std::endl;
	std::cerr << "ExtendedInt64q::max_numerator_bits=" <<
		CGAL::ExtendedInt64q<CGAL::Gmpq>::max_numerator_bits << std::endl;
	std::cerr << "ExtendedInt64q::max_denominator_bits=" <<
		CGAL::ExtendedInt64q<CGAL::Gmpq>::max_denominator_bits << std::endl;
#endif
}

template<>
TriangulationCreatorInExactIntersectionsConstrainedDelaunaySp::TriangulationCreatorConstrainedDelaunay(int significands, int intersectionSignificands) :
m_tr(
	TriangulationCreatorInExactIntersectionsConstrainedDelaunaySp::Tr::Geom_traits(
		LIB_RATSS_NAMESPACE::Conversion<TriangulationCreatorInExactIntersectionsConstrainedDelaunaySp::Tr::Geom_traits::FT>::moveFrom(
			LIB_RATSS_NAMESPACE::convert<mpq_class>(std::numeric_limits<uint32_t>::max()-1)/
			LIB_RATSS_NAMESPACE::convert<mpq_class>(std::numeric_limits<uint32_t>::max())
		),
		significands,
		intersectionSignificands
	)
)
{}

template<>
TriangulationCreatorInExactIntersectionsConstrainedDelaunaySpk::TriangulationCreatorConstrainedDelaunay(int significands, int intersectionSignificands) :
m_tr(
	TriangulationCreatorInExactIntersectionsConstrainedDelaunaySpk::Tr::Geom_traits(
		dts2::Kernel_sp::AuxiliaryPointsGenerator(),
		significands,
		intersectionSignificands
	)
)
{}

template<>
TriangulationCreatorInExactIntersectionsConstrainedDelaunaySpk64::TriangulationCreatorConstrainedDelaunay(int significands, int intersectionSignificands) :
m_tr(
	TriangulationCreatorInExactIntersectionsConstrainedDelaunaySpk64::Tr::Geom_traits(
		dts2::Kernel_sp_64::AuxiliaryPointsGenerator(),
		significands,
		intersectionSignificands
	)
)
{}

class Config: public ratss::BasicCmdLineOptions {
public:
	TriangulationType triangType;
	GraphOutputType got = GOT_INVALID;
	GraphInputType git = GIT_INVALID;
	TriangulationInputOrder tio = TIO_NODES_EDGES;
	int intersectSignificands = -1;
	std::string rewriteFileName;
	bool check{false};
public:
	bool autoChangedOptions{false};
public:
	using ratss::BasicCmdLineOptions::parse;
public:
	virtual bool parse(const std::string & currentToken,int & i, int argc, char ** argv) override;
	virtual void parse_completed() override;
	void help(std::ostream & out) const;
	void print(std::ostream & out) const;
};

struct Data {
	using Points = std::vector<std::pair<Point3, VertexInfo>>;

	Data();
	~Data();

	ratss::FloatPoint ip;
	ratss::RationalPoint op;
	ratss::ProjectS2 proj;

	std::vector<std::pair<Point3, VertexInfo>> points;
	std::vector<std::pair<int, int>> edges;
	TriangulationCreator * tc;
	
	void init(const Config& cfg);
	void read(InputOutput& io, const Config& cfg);
	void rewrite(InputOutput&io, const Config& cfg);
	void create(InputOutput& io, const Config& cfg);
	void write(InputOutput& io, const Config& cfg);
	
	void info(InputOutput& io);
	
	void readNodesEdges(InputOutput& io, const Config& cfg);
	void readEdges(InputOutput& io, const Config& cfg);
	Point3 readPoint(std::istream& is, const Config& cfg);
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
	
	if (cfg.autoChangedOptions) {
		io.info() << "NOTICE: Some options were added/changed to accomodate the requested type of triangulation." << std::endl;
	}

	if (cfg.verbose) {
		cfg.print(io.info());
		io.info() << std::endl;
	}
	TimeMeasurer ttm;
	TimeMeasurer tm;
	MemUsage mem;
	
	ttm.begin();
	
	io.info() << "Initializing triangulation..." << std::endl;
	data.init(cfg);
	
	{
		io.info() << "Reading graph..." << std::endl;
		tm.begin();
		data.read(io, cfg);
		tm.end();
		::malloc_trim(0);
		mem.update();
		io.info() << "Memory usage: " << mem << std::endl;
		io.info() << "Reading the graph took " << tm << std::endl;
	}
	
	if (cfg.git != GIT_NODES_EDGES_BY_POINTS_BINARY && cfg.tio == TIO_NODES_EDGES_BY_POINTS_BINARY) {
		io.info() << "Rewriting the graph..." << std::endl;
		tm.begin();
		data.rewrite(io, cfg);
		tm.end();
		io.info() << "Memory usage: " << mem << std::endl;
		io.info() << "Rewriting the graph took " << tm << std::endl;
		data.points = decltype(data.points)();
		data.edges = decltype(data.edges)();
	}
	
	{
		io.info() << "Creating triangulation..." << std::endl;
		tm.begin();
		data.create(io, cfg);
		tm.end();
		::malloc_trim(0);
		mem.update();
		data.info(io);
		io.info() << "Memory usage: " << mem << std::endl;
		io.info() << "Time: " << tm << std::endl;
	}
	
	
	io.info() << "Writing triangulation..." << std::endl;
	data.write(io, cfg);
	
	ttm.end();
	io.info() << "Total time: " << ttm << std::endl;
	
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
		else if (type == "d64" || type == "delaunay-64") {
			triangType = TT_DELAUNAY_64;
		}
		else if (type == "ds" || type == "delaunay-spherical") {
			triangType = TT_DELAUNAY_SPHERICAL;
		}
		else if (type == "chi" || type == "convexhull-inexact") {
			triangType = TT_CONVEX_HULL_INEXACT;
		}
		else if (type == "ch" || type == "convexhull") {
			triangType = TT_CONVEX_HULL;
		}
		else if (type == "ch64" || type == "convexhull-64") {
			triangType = TT_CONVEX_HULL_64;
		}
		else if (type == "c" || type == "constrained") {
			triangType = TT_CONSTRAINED;
		}
		else if (type == "cx" || type == "constrained-intersection") {
			triangType = TT_CONSTRAINED_INEXACT;
		}
		else if (type == "cx64" || type == "constrained-intersection-64") {
			triangType = TT_CONSTRAINED_INEXACT_64;
		}
		else if (type == "cxsp" || type == "constrained-intersection-sp") {
			triangType = TT_CONSTRAINED_INEXACT_SP;
		}
		else if (type == "cxspk" || type == "constrained-intersection-spk") {
			triangType = TT_CONSTRAINED_INEXACT_SPK;
		}
		else if (type == "cxspk64" || type == "constrained-intersection-spk-64") {
			triangType = TT_CONSTRAINED_INEXACT_SPK64;
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
	else if (token == "-go" && i+1 < argc) {
		std::string type( argv[i+1] );
		if (type == "none") {
			got = GOT_NONE;
		}
		else if (type == "wx") {
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
	else if (token == "-gi" && i+1 < argc) {
		std::string type( argv[i+1] );
		if (type == "ne" || type == "nodes-edges") {
			git = GIT_NODES_EDGES;
		}
		else if (type == "e" || type == "edges") {
			git = GIT_EDGES;
		}
		else {
			throw ratss::BasicCmdLineOptions::ParseError("Unknown graph input type: " + type);
		}
		++i;
	}
	else if (token == "-io" && i+1 < argc) {
		std::string type( argv[i+1] );
		if (type == "ne" || type == "nodes-edges") {
			tio= TIO_NODES_EDGES;
		}
		else if (type == "e" || type == "edges") {
			tio = TIO_EDGES;
		}
		else if (type == "nei" || type == "nodes-edges-iterative") {
			tio = TIO_NODES_EDGES_BY_POINTS_BINARY;
		}
		else {
			throw ratss::BasicCmdLineOptions::ParseError("Unknown graph input type: " + type);
		}
		++i;
	}
	else if (token == "-is" && i+1 < argc) {
		intersectSignificands = ::atoi(argv[i+1]);
		++i;
	}
	else if (token == "--check") {
		check = true;
	}
	else {
		return false;
	}
	return true;
}

void Config::parse_completed() {
	if (intersectSignificands < 2) {
		intersectSignificands = significands;
		autoChangedOptions = true;
	}
	if (triangType == TT_CONSTRAINED_INEXACT_64 || triangType == TT_CONSTRAINED_INEXACT_SP || triangType == TT_CONSTRAINED_INEXACT_SPK || triangType == TT_CONSTRAINED_INEXACT_SPK64 || triangType == TT_CONVEX_HULL_64 || triangType == TT_DELAUNAY_64) {
		significands = 31;
		intersectSignificands = 31;
		snapType = ratss::ST_PLANE | ratss::ST_FX | ratss::ST_NORMALIZE;
		autoChangedOptions = true;
		
	}
	if (triangType == TT_CONVEX_HULL_64 || triangType == TT_CONVEX_HULL || triangType == TT_CONVEX_HULL_INEXACT) {
		if (tio != TIO_NODES_EDGES && tio != TIO_NODES_EDGES_BY_POINTS_BINARY) {
			tio = TIO_NODES_EDGES;
			autoChangedOptions = true;
		}
	}
	if (tio == TIO_NODES_EDGES_BY_POINTS_BINARY && rewriteFileName.empty()) {
		rewriteFileName = inFileName + ".triang.temp";
	}
}

void Config::help(std::ostream & out) const {
	out << "triang OPTIONS:\n"
		"\t-t type\ttype = [d,delaunay,d64,delaunay-64,ds,delaunay-spherical,chi,convexhull-inexact,ch,convexhull,ch64,convexhull-64,c,constrained,cx,constrained-intersection,cx64,constrained-intersection-64,cxsp,constrained-intersection-sp,cxspk-64,constrained-intersection-spk-64,cxe,constrained-intesection-exact, cxs, constrained-intersection-exact-spherical]\n"
		"\t-go type\tgraph output type = [none, wx, witout_special, simplest, simplest_andre]\n"
		"\t-gi type\tgraph input type = [ne, nodes-edges, e, edges, nei, nodes-edges-iterative]\n"
		"\t-io type\tinput order type = [ne, nodes-edges, e, edges, nei, nodes-edges-iterative]\n"
		"\t-is num\tsignificands used to calculate intersection points\n"
		"\t--check\tCheck points and triangulation\n";
	ratss::BasicCmdLineOptions::options_help(out);
	out << '\n';
	out << "For gi=ne the input format is as follows:\n"
		"\tnumber of nodes\n"
		"\tnumber of edges\n"
		"\tnodes - one node per line\n"
		"\tedges - one edge per line\n"
		"\tA node is any coordinate type supported by the -if switch\n"
		"\tThe i-th node gets the id i. There shall not be multiple nodes with the same coordinates.\n"
		"\tAn edge has the format: source-node-id target-node-id\n"
		"For gi=e the input format is:\n"
		"\tsrc-node target-node\n"
		"For gi=nei the input format is:\n"
		"\tcomplicated binary format, code IS documentation!";
	out << std::endl;
}

void Config::print(std::ostream & out) const {
	out << "Triangulation type: ";
	switch (triangType) {
	case TT_DELAUNAY:
		out << "delaunay";
		break;
	case TT_DELAUNAY_64:
		out << "delaunay using ExtendedInt64 kernel";
		break;
	case TT_DELAUNAY_SPHERICAL:
		out << "delaunay using spherical kernel";
		break;
	case TT_CONVEX_HULL:
		out << "convex hull";
		break;
	case TT_CONVEX_HULL_INEXACT:
		out << "convex hull inexact";
		break;
	case TT_CONVEX_HULL_64:
		out << "convex hull using ExtendedInt64 kernel";
		break;
	case TT_CONSTRAINED:
		out << "constrained no intersections";
		break;
	case TT_CONSTRAINED_INEXACT:
		out << "constrained in-exact intersections";
		break;
	case TT_CONSTRAINED_INEXACT_64:
		out << "constrained in-exact intersections using ExtendedInt64 kernel";
		break;
	case TT_CONSTRAINED_INEXACT_SP:
		out << "constrained in-exact intersections using ExtendedInt64 kernel and Point_sp";
		break;
	case TT_CONSTRAINED_INEXACT_SPK:
		out << "constrained in-exact intersections using Kernel_sp";
		break;
	case TT_CONSTRAINED_INEXACT_SPK64:
		out << "constrained in-exact intersections using Kernel_sp with ExtendedInt64";
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
	out << "Graph input type: ";
	switch (git) {
	case GIT_NODES_EDGES:
		out << "nodes and edges";
		break;
	case GIT_EDGES:
		out << "edges";
		break;
	default:
		out << "invalid";
		break;
	};
	out << '\n';
	out << "Graph output type: ";
	switch (got) {
	case GOT_NONE:
		out << "none";
		break;
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
	out << "Input order: ";
	switch (tio) {
	case TIO_NODES_EDGES:
		out << "nodes then edges";
		break;
	case TIO_EDGES:
		out << "edges";
		break;
	case TIO_NODES_EDGES_BY_POINTS_BINARY:
		out << "nodes then edges iterative";
		break;
	default:
		out << "invalid";
		break;
	};
	out << '\n';
	if (triangType == TT_CONSTRAINED_INEXACT || triangType == TT_CONSTRAINED_INEXACT_64 || triangType == TT_CONSTRAINED_INEXACT_SP || triangType == TT_CONSTRAINED_INEXACT_SPK64 || triangType == TT_CONSTRAINED_INEXACT_SPK64) {
		out << "Intersection point significands: " << intersectSignificands << '\n';
	}
	std::cout << "Check: " << (check ? "yes" : "no") << '\n';
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
		tc = new TriangulationCreatorDelaunayEpeck(cfg.significands);
		break;
	case TT_DELAUNAY_64:
		tc = new TriangulationCreatorDelaunayFsceik(cfg.significands);
		break;
	case TT_DELAUNAY_SPHERICAL:
		tc = new TriangulationCreatorDelaunaySpherical(cfg.significands);
		break;
	case TT_CONVEX_HULL_INEXACT:
		tc = new EpickConvexHullTriangulationCreator();
		break;
	case TT_CONVEX_HULL:
		tc = new EpeckConvexHullTriangulationCreator();
		break;
	case TT_CONVEX_HULL_64:
		tc = new Ei64ConvexHullTriangulationCreator();
		break;
	case TT_CONSTRAINED:
		tc = new TriangulationCreatorNoIntersectionsConstrainedDelaunay(cfg.significands);
		break;
	case TT_CONSTRAINED_INEXACT:
		tc = new TriangulationCreatorInExactIntersectionsConstrainedDelaunay(cfg.significands, cfg.intersectSignificands);
		break;
	case TT_CONSTRAINED_INEXACT_64:
		tc = new TriangulationCreatorInExactIntersectionsConstrainedDelaunay64(cfg.significands, cfg.intersectSignificands);
		break;
	case TT_CONSTRAINED_INEXACT_SP:
		tc = new TriangulationCreatorInExactIntersectionsConstrainedDelaunaySp(cfg.significands, cfg.intersectSignificands);
		break;
	case TT_CONSTRAINED_INEXACT_SPK:
		tc = new TriangulationCreatorInExactIntersectionsConstrainedDelaunaySpk(cfg.significands, cfg.intersectSignificands);
		break;
	case TT_CONSTRAINED_INEXACT_SPK64:
		tc = new TriangulationCreatorInExactIntersectionsConstrainedDelaunaySpk64(cfg.significands, cfg.intersectSignificands);
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
	switch (cfg.git) {
	case GIT_NODES_EDGES:
		readNodesEdges(io, cfg);
		break;
	case GIT_EDGES:
		readEdges(io, cfg);
	case GIT_NODES_EDGES_BY_POINTS_BINARY:
		break;
	default:
		throw std::runtime_error("Invalid graph input format");
		break;
	}
}


class MyPropertyMap {
public:
	using BaseTrait = Point3::K;
	using key_type = uint32_t;
	using value_type = BaseTrait::Point_3;
	using reference = value_type;
	using category = boost::readable_property_map_tag;
public:
	MyPropertyMap(Data::Points * pts) : m_pts(pts) {}
	value_type get(key_type p) const {
		return (value_type) m_pts->at(p).first;
	}
private:
	Data::Points * m_pts;
};

MyPropertyMap::value_type get(const MyPropertyMap & p, MyPropertyMap::key_type k) {
	return p.get(k);
}

MyPropertyMap::value_type get(MyPropertyMap * p, MyPropertyMap::key_type k) {
	return p->get(k);
}

void Data::rewrite(InputOutput& io, const Config& cfg) {
	assert(cfg.tio == TIO_NODES_EDGES_BY_POINTS_BINARY);

	using SortTraits = CGAL::Spatial_sort_traits_adapter_3<MyPropertyMap::BaseTrait, MyPropertyMap>;
	
	InputOutput myio;
	myio.setOutput(cfg.rewriteFileName, std::ios_base::out | std::ios_base::binary);
	
	BinaryIo bio;
	bio.setIo(&myio);
	
	std::vector<MyPropertyMap::key_type> remap(points.size(), 0); //maps new order to old order
	for(MyPropertyMap::key_type i(0), s(points.size()); i < s; ++i) {
		remap.at(i) = i;
	}
	
	CGAL::spatial_sort(remap.begin(), remap.end(), SortTraits( MyPropertyMap(&points) ));
	
	//now write the points and edges to storage
	bio.write<uint32_t>(points.size());
	bio.write<uint32_t>(edges.size());
	for(uint32_t i(0), s(remap.size()); i < s; ++i) {
		bio.write<Point3>(points.at(remap.at(i)).first);
	}
	
	//we now need the inverse of remap
	std::vector<uint32_t> remap_inverse(remap.size(), 0); //maps old order to new order
	for(uint32_t i(0), s(remap.size()); i < s; ++i) {
		remap_inverse.at(remap.at(i)) = i;
	}
	
	std::sort(edges.begin(), edges.end(),
		[&remap_inverse](const std::pair<int, int> & a, const std::pair<int, int> & b) {
			return remap_inverse.at(a.first) < remap_inverse.at(b.first);
		}
	);
	
	for(const auto & e : edges) {
		bio.write<Point3>(points.at(e.first).first);
		bio.write<Point3>(points.at(e.second).first);
	}
	io.output().flush();
}


void Data::readNodesEdges(InputOutput & io, const Config & cfg) {
	std::istream & is = io.input();

	std::size_t num_points, num_edges;
	is >> num_points >> num_edges;
	
	io.info() << "Need to fetch " << num_points << " points and " << num_edges << " edges" << std::endl;
	
	points.reserve(num_points);
	edges.reserve(num_edges);
	
	if (cfg.progress) {
		io.info() << std::endl;
	}
	std::size_t counter = 0;
	while( is.good() && points.size() < num_points) {
		if (is.peek() == '\n') {
			is.get();
			continue;
		}
		points.emplace_back(readPoint(is, cfg), VertexInfo(counter));
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

void Data::readEdges(InputOutput& io, const Config& cfg) {
	ratss::FloatPoint ip;
	ratss::RationalPoint op;
	
	std::istream & is = io.input();
	
	if (cfg.progress) {
		io.info() << std::endl;
	}
	std::size_t counter = 0;
	while( is.good() ) {
		if (is.peek() == '\n') {
			is.get();
			continue;
		}
		Point3 src( readPoint(is, cfg) );
		Point3 tgt( readPoint(is, cfg) );
		
		int srcId = points.size();
		int tgtId = srcId+1;
		
		points.emplace_back(src, VertexInfo(srcId));
		points.emplace_back(tgt, VertexInfo(tgtId));
		edges.emplace_back(srcId, tgtId);
		
		++counter;
		if (cfg.progress && counter % 1000 == 0) {
			io.info() << '\xd' << counter/1000 << "k" << std::flush;
		}
	}
	if (cfg.progress) {
		io.info() << std::endl;
	}
}

Point3 Data::readPoint(std::istream& is, const Config& cfg) {
	bool geoFormat = (cfg.inFormat & (ratss::FloatPoint::FM_GEO|ratss::FloatPoint::FM_SPHERICAL));
	ratss::PointBase::Format format = (geoFormat ? ratss::PointBase::FM_CARTESIAN_FLOAT128 : cfg.inFormat);
	bool opFromIp = !cfg.rationalPassThrough;
	if (cfg.rationalPassThrough) {
		op.assign(is, format, cfg.precision, 3);
		if (!op.valid()) {
			ip.assign(op.coords.begin(), op.coords.end(), cfg.precision);
			opFromIp = true;
		}
	}
	else {
		ip.assign(is, format, cfg.precision, 3);
	}
	if (opFromIp) {
		op.clear();
		op.resize(3);
		ip.setPrecision(cfg.precision);
		switch(cfg.inFormat) {
			case ratss::FloatPoint::FM_GEO:
				if (ip.coords.size() != 2) {
					throw std::runtime_error("Invalid input point");
				}
				proj.projectFromGeo(ip.coords[0], ip.coords[1], op.coords[0], op.coords[1], op.coords[2], cfg.significands, cfg.snapType);
				break;
			case ratss::FloatPoint::FM_SPHERICAL:
				if (ip.coords.size() != 2) {
					throw std::runtime_error("Invalid input point");
				}
				proj.projectFromSpherical(ip.coords[0], ip.coords[1], op.coords[0], op.coords[1], op.coords[2], cfg.significands, cfg.snapType);
				break;
			default:
				proj.snap(ip.coords.begin(), ip.coords.end(), op.coords.begin(), cfg.snapType, cfg.significands);
				break;
		}
	}
	assert(op.valid());
	
	if (cfg.check) {
		if (!op.valid()) {
			throw std::runtime_error("Invalid point detected");
		}
	}
	
	using FT = Point3::FT;
	FT x = ratss::Conversion<FT>::moveFrom( std::move(op.coords.at(0)) );
	FT y = ratss::Conversion<FT>::moveFrom( std::move(op.coords.at(1)) );
	FT z = ratss::Conversion<FT>::moveFrom( std::move(op.coords.at(2)) );
	
	Point3 p3(x, y, z);
	
	return p3;
}

void Data::create(InputOutput& io, const Config& cfg) {
	if (cfg.tio == TIO_NODES_EDGES) {
		tc->create(points, edges, io, true);
	}
	else if (cfg.tio == TIO_EDGES) {
		if (cfg.progress) {
			io.info() << '\n' << std::flush;
		}
		
		for(
			std::size_t i(0),
			s(edges.size())
			;i < s; ++i)
		{
			const auto & e = edges[i];
			const auto & src = points.at(e.first).first;
			const auto & tgt = points.at(e.second).first;
			tc->add(src, tgt);
			if (cfg.progress && i % 1000 == 0) {
				io.info() << '\xd' << i/1000 << 'k';
			}
		}
		
		if (cfg.progress) {
			io.info() << std::endl;
		}
	}
	else if (cfg.tio == TIO_NODES_EDGES_BY_POINTS_BINARY) {
		InputOutput myio;
		BinaryIo bio;
		
		if (cfg.git == GIT_NODES_EDGES_BY_POINTS_BINARY) {
			io.input();
			assert(false);
		}
		else {
			myio.setInput(cfg.rewriteFileName, std::ios_base::in | std::ios_base::binary);
			bio.setIo(&myio);
		}
		
		uint32_t nc = bio.get<uint32_t>();
		uint32_t ec = bio.get<uint32_t>();
		
		if (cfg.tio == TIO_NODES_EDGES_BY_POINTS_BINARY) {
			io.info() << "Inserting points..." << std::flush;
			TimeMeasurer tm;
			tm.begin();
			for(uint32_t i(0); i < nc; ++i) {
				Point3 p = bio.get<Point3>();
				tc->add(p);
			}
			tm.end();
			io.info() << tm << std::endl;
			io.info() << "Inserting constraints..." << std::flush;
			tm.begin();
			for(uint32_t i(0); i < ec; ++i) {
				auto p1 = bio.get<Point3>();
				auto p2 = bio.get<Point3>();
				tc->add(p1, p2);
			}
			io.info() << tm << std::endl;
		}
	}
	points = decltype(points)();
	edges = decltype(edges)();
}

void Data::info(InputOutput& io) {
	tc->info(io);
}

void Data::write(InputOutput & io, const Config & cfg) {
	tc->got = cfg.got;
	tc->pointFormat = cfg.outFormat;
	tc->write(io);
}
