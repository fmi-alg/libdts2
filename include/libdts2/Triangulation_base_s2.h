#ifndef LIB_DTS2_TRIANGULATION_BASE_S2_H
#define LIB_DTS2_TRIANGULATION_BASE_S2_H
#pragma once

#include <libdts2/constants.h>
#include <libratss/GeoCoord.h>
#include <libratss/SphericalCoord.h>
#include <libratss/Conversion.h>
#include <libratss/enum.h>

#include <CGAL/internal/info_check.h>
#include <CGAL/enum.h>

#include <functional>
#include <type_traits>
#include <iterator>
#include <assert.h>

//BUG: check if z-coordinate is out of range after projection!

namespace LIB_DTS2_NAMESPACE {
namespace internal {
namespace Triangulation_base_s2 {
//based on https://rmf.io/cxx11/almost-static-if/
	enum class enabler {};
}}

///This is the base class for all triangulations on the sphere
template<typename T_BASE_TRIANGULATION>
class Triangulation_base_s2: private T_BASE_TRIANGULATION {
public:
	using MyBaseClass = T_BASE_TRIANGULATION;
public:
	using Triangulation = T_BASE_TRIANGULATION;
	using Triangulation_data_structure = typename Triangulation::Triangulation_data_structure;
	using Geom_traits = typename Triangulation::Geom_traits;
	using Trs = Triangulation;
	using Tds = Triangulation_data_structure;
public:
	using typename Triangulation::Locate_type;
	using Triangulation::VERTEX;
	using Triangulation::EDGE;
	using Triangulation::FACE;
	using Triangulation::OUTSIDE_CONVEX_HULL;
	using Triangulation::OUTSIDE_AFFINE_HULL;
public:
	using FT = typename Geom_traits::FT;
	using Point = typename Geom_traits::Point;
	using Point_3 = typename Geom_traits::Point_3;
	using GeoCoord = LIB_RATSS_NAMESPACE::GeoCoord;
	using SphericalCoord = LIB_RATSS_NAMESPACE::SphericalCoord;
public: //triangulation types
	using Segment = typename Triangulation::Segment;
	using Edge = typename Triangulation::Edge;
	
	using Vertex = typename Triangulation::Vertex;
	using Vertex_handle = typename Triangulation::Vertex_handle;
	using Vertex_circulator = typename Triangulation::Vertex_circulator;

	using Face = typename Triangulation::Face;
	using Face_handle = typename Triangulation::Face_handle;
	using Face_circulator = typename Triangulation::Face_circulator;
	
	using Edge_circulator = typename Triangulation::Edge_circulator;
	
	using Finite_vertices_iterator = typename Triangulation::Finite_vertices_iterator;
	using Finite_edges_iterator = typename Triangulation::Finite_edges_iterator;
	using Finite_faces_iterator = typename Triangulation::Finite_faces_iterator;
public:
	using Project_on_sphere = typename Geom_traits::Project_on_sphere;
public:
	Triangulation_base_s2(int significands);
	Triangulation_base_s2(const Geom_traits & traits);
	Triangulation_base_s2(Triangulation_base_s2 && other);
	Triangulation_base_s2(const Triangulation_base_s2& other) = delete;
	Triangulation_base_s2 & operator=(Triangulation_base_s2 && other);
public: //insertion
	Vertex_handle insert(const Point & p, const Face_handle & fh = Face_handle());
	Vertex_handle insert(double lat, double lon, const Face_handle & fh = Face_handle());
	Vertex_handle insert(const std::pair<double, double> & p, const Face_handle & fh = Face_handle());
	Vertex_handle insert(const SphericalCoord & p, const Face_handle & fh = Face_handle());
	Vertex_handle insert(const GeoCoord & p, const Face_handle & fh = Face_handle());

	template<typename T_ITERATOR>
	void insert(T_ITERATOR begin, T_ITERATOR end) { insert(begin, end, true); }

	template<typename T_ITERATOR>
	void insert(T_ITERATOR begin, T_ITERATOR end, bool snap);
	
	template<typename T_ITERATOR>
	void insert_with_info(T_ITERATOR begin, T_ITERATOR end) { insert_with_info(begin, end, true); }
	
	template<typename T_ITERATOR>
	void insert_with_info(T_ITERATOR begin, T_ITERATOR end, bool snap);
public:
	using MyBaseClass::locate;
public: //removal
	void remove(const Vertex_handle & vh);
public:
	Finite_vertices_iterator finite_vertices_begin() const;
	Finite_vertices_iterator finite_vertices_end() const;

	Finite_edges_iterator finite_edges_begin() const;
	Finite_edges_iterator finite_edges_end() const;
	
	Finite_faces_iterator finite_faces_begin() const;
	Finite_faces_iterator finite_faces_end() const;
	
	Vertex_circulator incident_vertices(const Vertex_handle & vh) const;
	Edge_circulator incident_edges(const Vertex_handle & vh) const;
	Face_circulator incident_faces(const Vertex_handle & vh) const;
public:
	Segment segment(const Edge & e) const;
public: //size info
	std::size_t number_of_vertices() const;
	std::size_t number_of_edges();
	std::size_t number_of_faces() const;
public: //information stuff
	bool is_infinite(const Vertex_handle & vh) const;
	bool is_auxiliary(const Vertex_handle & vh) const;
	/// the same as is_infinite || is_auxiliary
	bool is_special(const Vertex_handle & vh) const;
	std::size_t degree(const Vertex_handle & vh) const;
	bool is_valid();
public: //helpers
	int ccw(int v) const;
	int cw(int v) const;
public: //projection
	Point_3 project(const std::pair<double, double> & latLon) const;
	Point_3 project(double lat, double lon) const;
	Point_3 project(const GeoCoord & geoCoord) const;
	Point_3 project(const SphericalCoord & sphericalCoord) const;
	Point_3 project(const Point_3 & p) const;
	Point_3 projectFromSpherical(const mpfr::mpreal & theta, const mpfr::mpreal & phi) const;
	SphericalCoord toSpherical(const Point_3 & p) const;
	GeoCoord toGeo(const Point_3 & p) const;
public:
	int significands() const;
	const Geom_traits & geom_traits() const;
	const Triangulation & trs() const;
	const Tds & tds() const;
	//the z-coordinate of points that are invalid since they are within the aux-triangle
	const FT & epsZ() const;
protected:

	template<typename T_ITERATOR, typename T_ITERATOR_VALUE_TYPE, bool T_SNAP>
	void insert_with_info_impl(T_ITERATOR begin, T_ITERATOR end, T_ITERATOR_VALUE_TYPE* /*dummy*/);
	
	template<typename T_ITERATOR, typename T_ITERATOR_VALUE_TYPE, bool T_SNAP>
	void insert_with_info_impl(T_ITERATOR begin, T_ITERATOR end,
		typename std::enable_if<
			std::is_same<
				T_ITERATOR_VALUE_TYPE,
				std::pair<Point_3, typename CGAL::internal::Info_check<Vertex>::type>
			>::value,
			T_ITERATOR_VALUE_TYPE
		>::type * /*dummy*/);
	
	template<typename T_ITERATOR, typename T_ITERATOR_VALUE_TYPE, bool T_SNAP>
	void insert_impl(T_ITERATOR begin, T_ITERATOR end, T_ITERATOR_VALUE_TYPE* /*dummy*/);
	
	template<typename T_ITERATOR, typename T_ITERATOR_VALUE_TYPE, bool T_SNAP>
	void insert_impl(T_ITERATOR begin, T_ITERATOR end,
		typename std::enable_if<
			std::is_same<
				T_ITERATOR_VALUE_TYPE,
				Point_3
			>::value,
			T_ITERATOR_VALUE_TYPE
		>::type * /*dummy*/);
	
	template<typename T_ITERATOR, typename T_ITERATOR_VALUE_TYPE, bool T_SNAP>
	void insert_impl(T_ITERATOR begin, T_ITERATOR end,
		typename std::enable_if<
			std::is_same<
				T_ITERATOR_VALUE_TYPE,
				std::pair<Point_3, typename CGAL::internal::Info_check<Vertex>::type>
			>::value,
			T_ITERATOR_VALUE_TYPE
		>::type * /*dummy*/);
	
	template<typename T_ITERATOR, typename T_ITERATOR_VALUE_TYPE, bool T_SNAP>
	void insert_impl(T_ITERATOR begin, T_ITERATOR end,
		typename std::enable_if<
			std::is_same<
				T_ITERATOR_VALUE_TYPE,
				std::pair<SphericalCoord, typename CGAL::internal::Info_check<Vertex>::type>
			>::value,
			T_ITERATOR_VALUE_TYPE
		>::type * /*dummy*/);
	
	template<typename T_ITERATOR, typename T_ITERATOR_VALUE_TYPE, bool T_SNAP>
	void insert_impl(T_ITERATOR begin, T_ITERATOR end,
		typename std::enable_if<
			std::is_same<
				T_ITERATOR_VALUE_TYPE,
				std::pair<GeoCoord, typename CGAL::internal::Info_check<Vertex>::type>
			>::value,
			T_ITERATOR_VALUE_TYPE
		>::type * /*dummy*/);
protected:
	void selfCheck();
protected:
	Triangulation & trs();
	Tds & tds();
	const Project_on_sphere & proj() const;
private:
	void addAuxiliaryPoints();
	Point_3 generateAuxPoint(FT xp, FT yp) const;
private:
	Triangulation m_cdts;
	Project_on_sphere m_p;
};

} //end namespace LIB_DTS2_NAMESPACE

//and now the implementation
namespace LIB_DTS2_NAMESPACE {

#define TMPL_HDR template<typename T_BASE_TRIANGULATION>
#define TMPL_CLS Triangulation_base_s2<T_BASE_TRIANGULATION>

//BEGIN constructors
TMPL_HDR
TMPL_CLS::Triangulation_base_s2(int _significands) :
Triangulation_base_s2(
	Geom_traits(
		LIB_RATSS_NAMESPACE::Conversion<FT>::moveFrom(mpq_class(std::numeric_limits<uint64_t>::max()-1, std::numeric_limits<uint64_t>::max())),
		_significands
	)
)
{}

TMPL_HDR
TMPL_CLS::Triangulation_base_s2(const Geom_traits & traits) :
MyBaseClass(traits),
m_cdts(*this),
m_p( this->geom_traits().project_on_sphere_object() )
{
	if (! (epsZ() > 0 && epsZ() < 1) ) {
		throw std::runtime_error("libdts2::Triangulation_base_s2: epsZ musst be smaller than 1 and larger than 0");
	}
	addAuxiliaryPoints();
	selfCheck();
}

TMPL_HDR
TMPL_CLS::Triangulation_base_s2(Triangulation_base_s2 && other) :
MyBaseClass(other),
m_cdts(*this),
m_p(std::move(other.m_p))
{}

TMPL_HDR
TMPL_CLS & 
TMPL_CLS::operator=(TMPL_CLS && other) {
	MyBaseClass::operator=( std::move(other.m_cdts) );
	m_cdts = *this;
	m_p = std::move(other.m_p);
	return *this;
}
//END constructors

//BEGIN insertion operations

TMPL_HDR
typename TMPL_CLS::Vertex_handle
TMPL_CLS::insert(const Point & p, const Face_handle & fh) {
	if (p.x() * p.x() + p.y()*p.y() + p.z() * p.z() == 1) {
		return m_cdts.insert(p, fh);
	}
	else {
		return m_cdts.insert(project(p), fh);
	}
}

TMPL_HDR
typename TMPL_CLS::Vertex_handle 
TMPL_CLS::insert(double lat, double lon, const Face_handle & fh) {
	return insert(GeoCoord(lat, lon), fh);
}

TMPL_HDR
typename TMPL_CLS::Vertex_handle 
TMPL_CLS::insert(const std::pair<double, double> & latLon, const Face_handle & fh) {
	return insert(GeoCoord(latLon.first, latLon.second), fh);
}

TMPL_HDR
typename TMPL_CLS::Vertex_handle 
TMPL_CLS::insert(const SphericalCoord & p, const Face_handle & fh) {
	return m_cdts.insert(project(p), fh);
}

TMPL_HDR
typename TMPL_CLS::Vertex_handle 
TMPL_CLS::insert(const GeoCoord & p, const Face_handle & fh) {
	return m_cdts.insert(project(p), fh);
}

TMPL_HDR
template<typename T_ITERATOR>
void 
TMPL_CLS::insert(T_ITERATOR begin, T_ITERATOR end, bool snap) {
	using iterator_value_type = typename std::iterator_traits<T_ITERATOR>::value_type;
	if (snap) {
		insert_impl<T_ITERATOR, iterator_value_type, true>(begin, end, static_cast<iterator_value_type*>(0));
	}
	else {
		insert_impl<T_ITERATOR, iterator_value_type, false>(begin, end, static_cast<iterator_value_type*>(0));
	}
}

TMPL_HDR
template<typename T_ITERATOR>
void 
TMPL_CLS::insert_with_info(T_ITERATOR begin, T_ITERATOR end, bool snap) {
	using iterator_value_type = typename std::iterator_traits<T_ITERATOR>::value_type;
	if (snap) {
		insert_with_info_impl<T_ITERATOR, iterator_value_type, true>(begin, end, static_cast<iterator_value_type*>(0));
	}
	else {
		insert_with_info_impl<T_ITERATOR, iterator_value_type, false>(begin, end, static_cast<iterator_value_type*>(0));
	}
}
//END insertion operations

//BEGIN removal operations
TMPL_HDR
void
TMPL_CLS::remove(const Vertex_handle & vh) {
	m_cdts.remove(vh);
}
//END removal operations

//BEGIN iterators

TMPL_HDR
typename TMPL_CLS::Finite_vertices_iterator 
TMPL_CLS::finite_vertices_begin() const {
	return m_cdts.finite_vertices_begin();
}

TMPL_HDR
typename TMPL_CLS::Finite_vertices_iterator 
TMPL_CLS::finite_vertices_end() const {
	return m_cdts.finite_vertices_end();
}

TMPL_HDR
typename TMPL_CLS::Finite_edges_iterator 
TMPL_CLS::finite_edges_begin() const {
	return m_cdts.finite_edges_begin();
}

TMPL_HDR
typename TMPL_CLS::Finite_edges_iterator 
TMPL_CLS::finite_edges_end() const {
	return m_cdts.finite_edges_end();
}

TMPL_HDR
typename TMPL_CLS::Finite_faces_iterator 
TMPL_CLS::finite_faces_begin() const {
	return m_cdts.finite_faces_begin();
}

TMPL_HDR
typename TMPL_CLS::Finite_faces_iterator 
TMPL_CLS::finite_faces_end() const {
	return m_cdts.finite_faces_end();
}

TMPL_HDR
typename TMPL_CLS::Vertex_circulator 
TMPL_CLS::incident_vertices(const Vertex_handle & vh) const {
	return m_cdts.incident_vertices(vh);
}

TMPL_HDR
typename TMPL_CLS::Edge_circulator 
TMPL_CLS::incident_edges(const Vertex_handle & vh) const {
	return m_cdts.incident_edges(vh);
}

TMPL_HDR
typename TMPL_CLS::Face_circulator 
TMPL_CLS::incident_faces(const Vertex_handle & vh) const {
	return m_cdts.incident_faces(vh);
}
//END iterators

//BEGIN construction
TMPL_HDR
typename TMPL_CLS::Segment
TMPL_CLS::segment(const Edge & e) const {
	return m_cdts.segment(e);
}
//END construction

//BEGIN size info
TMPL_HDR
std::size_t 
TMPL_CLS::number_of_vertices() const {
	return m_cdts.number_of_vertices();
}

TMPL_HDR
std::size_t 
TMPL_CLS::number_of_edges() {
	return m_cdts.number_of_edges();
}

TMPL_HDR
std::size_t 
TMPL_CLS::number_of_faces() const {
	return m_cdts.number_of_faces();
}

//END size info
//BEGIN query functions
TMPL_HDR
bool 
TMPL_CLS::is_infinite(const Vertex_handle & vh) const {
	return m_cdts.is_infinite(vh);
}

TMPL_HDR
bool
TMPL_CLS::is_auxiliary(const Vertex_handle & vh) const {
	const auto & zc = vh->point().z();
	return zc == -1 || zc >= epsZ();
}

TMPL_HDR
bool
TMPL_CLS::is_special(const Vertex_handle & vh) const {
	return is_infinite(vh) || is_auxiliary(vh);
}

TMPL_HDR
std::size_t 
TMPL_CLS::degree(const Vertex_handle & vh) const {
	return m_cdts.degree(vh);
}

TMPL_HDR
bool
TMPL_CLS::is_valid() {
	return m_cdts.is_valid();
}
//END query functions
//BEGIN helpers

TMPL_HDR
int 
TMPL_CLS::ccw(int v) const {
	return m_cdts.ccw(v);
}

TMPL_HDR
int 
TMPL_CLS::cw(int v) const {
	return m_cdts.cw(v);
}

//END helpers
//BEGIN projection

TMPL_HDR
typename TMPL_CLS::Point_3 
TMPL_CLS::project(const std::pair<double, double> & latLon) const {
	return project(latLon.first, latLon.second);
}

TMPL_HDR
typename TMPL_CLS::Point_3 
TMPL_CLS::project(double lat, double lon) const {
	return project( GeoCoord(lat, lon) );
}

TMPL_HDR
typename TMPL_CLS::Point_3 
TMPL_CLS::project(const GeoCoord & geoCoord) const {
	return m_p(geoCoord);
}

TMPL_HDR
typename TMPL_CLS::Point_3 
TMPL_CLS::project(const SphericalCoord & sphericalCoord) const {
	return m_p(sphericalCoord);
}

TMPL_HDR
typename TMPL_CLS::Point_3
TMPL_CLS::project(const Point_3 & p) const {
	return m_p(p);
}

TMPL_HDR
typename TMPL_CLS::Point_3 
TMPL_CLS::projectFromSpherical(const mpfr::mpreal & /*theta*/, const mpfr::mpreal & /*phi*/) const {
	assert(false);
	return Point_3();
}

TMPL_HDR
LIB_RATSS_NAMESPACE::SphericalCoord 
TMPL_CLS::toSpherical(const Point_3 & p) const {
	SphericalCoord c;
	m_p.projector().toSpherical(p.x(), p.y(), p.z(), c.theta, c.phi, 64);
	return c;
}

TMPL_HDR
LIB_RATSS_NAMESPACE::GeoCoord 
TMPL_CLS::toGeo(const Point_3 & p) const {
	GeoCoord c;
	m_p.projector().toGeo(p.x(), p.y(), p.z(), c.lat, c.lon, 64);
	return c;
}
//END projection

//BEGIN member variable access
TMPL_HDR
int 
TMPL_CLS::significands() const {
	return geom_traits().significands();
}

TMPL_HDR
const typename TMPL_CLS::Geom_traits &
TMPL_CLS::geom_traits() const {
	return m_cdts.geom_traits();
}

TMPL_HDR
const typename TMPL_CLS::Triangulation & 
TMPL_CLS::trs() const {
	return m_cdts;
}

TMPL_HDR
const typename TMPL_CLS::Tds & 
TMPL_CLS::tds() const {
	return trs().tds();
}

TMPL_HDR
const typename TMPL_CLS::FT & 
TMPL_CLS::epsZ() const {
	return geom_traits().epsZ();
}

//END member variable access

//BEGIN private insertion functions

TMPL_HDR
template<typename T_ITERATOR, typename T_ITERATOR_VALUE_TYPE, bool T_SNAP>
void 
TMPL_CLS::insert_with_info_impl(T_ITERATOR begin, T_ITERATOR end, T_ITERATOR_VALUE_TYPE * /*dummy*/)
{
	std::vector< std::pair<Point_3, typename CGAL::internal::Info_check<Vertex>::type> > tmp;
	using std::distance;
	tmp.reserve(distance(begin, end));
	for(;begin != end; ++begin) {
		tmp.emplace_back(m_p(begin->first), begin->second);
	}
	m_cdts.insert(tmp.begin(), tmp.end());
}

TMPL_HDR
template<typename T_ITERATOR, typename T_ITERATOR_VALUE_TYPE, bool T_SNAP>
void 
TMPL_CLS::insert_with_info_impl(T_ITERATOR begin, T_ITERATOR end,
	typename std::enable_if<
		std::is_same<
			T_ITERATOR_VALUE_TYPE,
			std::pair<typename TMPL_CLS::Point_3, typename CGAL::internal::Info_check<Vertex>::type>
		>::value,
		T_ITERATOR_VALUE_TYPE
	>::type * /*dummy*/)
{
	if (T_SNAP) {
		std::vector< std::pair<Point_3, typename CGAL::internal::Info_check<Vertex>::type> > tmp;
		using std::distance;
		tmp.reserve(distance(begin, end));
		for(;begin != end; ++begin) {
			tmp.emplace_back(m_p(begin->first), begin->second);
		}
		m_cdts.insert(tmp.begin(), tmp.end());
	}
	else {
		m_cdts.insert(begin, end);
	}
}

TMPL_HDR
template<typename T_ITERATOR, typename T_ITERATOR_VALUE_TYPE, bool T_SNAP>
void 
TMPL_CLS::insert_impl(T_ITERATOR begin, T_ITERATOR end, T_ITERATOR_VALUE_TYPE * /*dummy*/)
{
	std::vector<Point_3> tmp;
	using std::distance;
	tmp.reserve( distance(begin, end) );
	for(;begin != end; ++begin) {
		tmp.emplace_back( m_p(*begin) );
	}
	m_cdts.insert(tmp.begin(), tmp.end());
}

TMPL_HDR
template<typename T_ITERATOR, typename T_ITERATOR_VALUE_TYPE, bool T_SNAP>
void 
TMPL_CLS::insert_impl(T_ITERATOR begin, T_ITERATOR end,
	typename std::enable_if<
		std::is_same<
			T_ITERATOR_VALUE_TYPE,
			Point_3
		>::value,
		T_ITERATOR_VALUE_TYPE
	>::type * /*dummy*/)
{
	if (T_SNAP) {
		std::vector<Point_3> tmp;
		using std::distance;
		tmp.reserve( distance(begin, end) );
		for(;begin != end; ++begin) {
			tmp.emplace_back( m_p(*begin) );
		}
		m_cdts.insert(tmp.begin(), tmp.end());
	}
	else {
		m_cdts.insert(begin, end);
	}
}

TMPL_HDR
template<typename T_ITERATOR, typename T_ITERATOR_VALUE_TYPE, bool T_SNAP>
void 
TMPL_CLS::insert_impl(T_ITERATOR begin, T_ITERATOR end,
	typename std::enable_if<
		std::is_same<
			T_ITERATOR_VALUE_TYPE,
			std::pair<typename TMPL_CLS::Point_3, typename CGAL::internal::Info_check<Vertex>::type>
		>::value,
		T_ITERATOR_VALUE_TYPE
	>::type * /*dummy*/)
{
	insert_with_info(begin, end, T_SNAP);
}

TMPL_HDR
template<typename T_ITERATOR, typename T_ITERATOR_VALUE_TYPE, bool T_SNAP>
void 
TMPL_CLS::insert_impl(T_ITERATOR begin, T_ITERATOR end,
	typename std::enable_if<
		std::is_same<
			T_ITERATOR_VALUE_TYPE,
			std::pair<SphericalCoord, typename CGAL::internal::Info_check<Vertex>::type>
		>::value,
		T_ITERATOR_VALUE_TYPE
	>::type * /*dummy*/)
{
	insert_with_info(begin, end, true);
}

TMPL_HDR
template<typename T_ITERATOR, typename T_ITERATOR_VALUE_TYPE, bool T_SNAP>
void 
TMPL_CLS::insert_impl(T_ITERATOR begin, T_ITERATOR end,
	typename std::enable_if<
		std::is_same<
			T_ITERATOR_VALUE_TYPE,
			std::pair<GeoCoord, typename CGAL::internal::Info_check<Vertex>::type>
		>::value,
		T_ITERATOR_VALUE_TYPE
	>::type * /*dummy*/)
{
	insert_with_info(begin, end, true);
}
//END private insertion functions

//BEGIN selfChecks
TMPL_HDR
void
TMPL_CLS::selfCheck() {
	GeoCoord p(40.0, 10.0); //lower left
	GeoCoord q(40.0, 12.0); //lower right
	GeoCoord r(42.0, 11.0); //upper mid
	GeoCoord s_in(41.0, 11.0);
	GeoCoord s_out(20.0, 11.0);
	
	Point_3 p3(project(p));
	Point_3 q3(project(q));
	Point_3 r3(project(r));
	Point_3 s3_in(project(s_in));
	Point_3 s3_out(project(s_out));

	{ //side of oriented circle tests

		
		auto soc = geom_traits().side_of_oriented_circle_2_object();
		
		assert(soc(p3, q3, r3, s3_in) == CGAL::LEFT_TURN);
		assert(soc(p3, q3, r3, s3_out) == CGAL::RIGHT_TURN);
	}
	{ // orientation2 test
		auto ot2 = geom_traits().orientation_2_object();
		
		assert(ot2(p3, q3, s3_in) == CGAL::LEFT_TURN);
		assert(ot2(q3, r3, s3_in) == CGAL::LEFT_TURN);
		assert(ot2(r3, p3, s3_in) == CGAL::LEFT_TURN);
		
// 		std::cout << geom_traits().doubleValue(p3.z()) << " " << geom_traits().doubleValue(q3.z()) << " " << geom_traits().doubleValue(r3.z()) << std::endl;
		
		assert(ot2(p3, q3, s3_out) == CGAL::RIGHT_TURN);
		assert(ot2(q3, r3, s3_out) == CGAL::LEFT_TURN);
	}
}
//END selfChecks

//BEGIN private member variable access

TMPL_HDR
typename TMPL_CLS::Triangulation & 
TMPL_CLS::trs() {
	return m_cdts;
}

TMPL_HDR
typename TMPL_CLS::Tds & 
TMPL_CLS::tds() {
	return trs().tds();
}

TMPL_HDR
const typename TMPL_CLS::Project_on_sphere &
TMPL_CLS::proj() const {
	return m_p;
}
//END private member variable access

//BEGIN private initial construction functions
TMPL_HDR
void 
TMPL_CLS::addAuxiliaryPoints() {
	
	//start value for generateAuxPoint
	FT epsilon(LIB_RATSS_NAMESPACE::Conversion<FT>::moveFrom(mpq_class(1, std::numeric_limits<uint16_t>::max())));
	auto p_1 = generateAuxPoint(0, epsilon);
	auto p_2 = generateAuxPoint(-epsilon, -epsilon/2);
	auto p_3 = generateAuxPoint(epsilon, -epsilon/2);
	
	if (p_1.z() < epsZ() || p_1.z() >= 1) {
		throw std::runtime_error("Construction auf auxiliary points failed. z-value=" + std::to_string(Geom_traits::doubleValue(p_1.z())));
	}
	if (p_2.z() < epsZ() || p_2.z() >= 1) {
		throw std::runtime_error("Construction auf auxiliary points failed. z-value=" + std::to_string(Geom_traits::doubleValue(p_2.z())));
	}
	if (p_3.z() < epsZ() || p_3.z() >= 1) {
		throw std::runtime_error("Construction auf auxiliary points failed. z-value=" + std::to_string(Geom_traits::doubleValue(p_3.z())));
	}
	
	tds().remove_first (trs().infinite_vertex());

	// Insert initial vertices and faces.

	auto vh_triangle_1 = tds().create_vertex ();
	auto vh_triangle_2 = tds().create_vertex ();
	auto vh_triangle_3 = tds().create_vertex ();
	auto vh_south_pole = tds().create_vertex ();
	auto vh_inf = tds().create_vertex ();

	auto fh_pyramid_13s = m_cdts.tds().create_face ();
	auto fh_pyramid_21s = m_cdts.tds().create_face ();
	auto fh_pyramid_32s = m_cdts.tds().create_face ();
	auto fh_inf_i31 = m_cdts.tds().create_face ();
	auto fh_inf_i12 = m_cdts.tds().create_face ();
	auto fh_inf_i23 = m_cdts.tds().create_face ();

	vh_triangle_1->set_point ( p_1 );
	vh_triangle_2->set_point ( p_2 );
	vh_triangle_3->set_point ( p_3 );
	
	vh_south_pole->set_point ( Point(0, 0, -1) );
	
	//put the infinte vertex outside of the triangulation (which is the north-pole)
	vh_inf->set_point ( Point(0, 0, 1) );
	
	vh_triangle_1->set_face (fh_pyramid_13s);
	vh_triangle_2->set_face (fh_pyramid_21s);
	vh_triangle_3->set_face (fh_pyramid_32s);
	vh_south_pole->set_face (fh_pyramid_13s);
	vh_inf->set_face (fh_inf_i31);

	fh_pyramid_13s->set_vertices (vh_triangle_1, vh_triangle_3, vh_south_pole);
	fh_pyramid_21s->set_vertices (vh_triangle_2, vh_triangle_1, vh_south_pole);
	fh_pyramid_32s->set_vertices (vh_triangle_3, vh_triangle_2, vh_south_pole);
	
	fh_inf_i31->set_vertices (vh_inf, vh_triangle_3, vh_triangle_1);
	fh_inf_i12->set_vertices (vh_inf, vh_triangle_1, vh_triangle_2);
	fh_inf_i23->set_vertices (vh_inf, vh_triangle_2, vh_triangle_3);
	
	fh_pyramid_13s->set_neighbors (fh_pyramid_32s, fh_pyramid_21s, fh_inf_i31);
	fh_pyramid_21s->set_neighbors (fh_pyramid_13s, fh_pyramid_32s, fh_inf_i12);
	fh_pyramid_32s->set_neighbors (fh_pyramid_21s, fh_pyramid_13s, fh_inf_i23);
	
	fh_inf_i31->set_neighbors (fh_pyramid_13s, fh_inf_i12, fh_inf_i23);
	fh_inf_i12->set_neighbors (fh_pyramid_21s, fh_inf_i23, fh_inf_i31);
	fh_inf_i23->set_neighbors (fh_pyramid_32s, fh_inf_i31, fh_inf_i12);

	m_cdts.set_infinite_vertex (vh_inf);
	m_cdts.tds().set_dimension (2);
	assert(m_cdts.tds().dimension() == 2);
	assert(m_cdts.tds().dimension() == m_cdts.dimension());
	assert(m_cdts.tds().dimension() == this->dimension());

	assert( is_valid() );
}

TMPL_HDR
typename TMPL_CLS::Point_3 
TMPL_CLS::generateAuxPoint(FT xp, FT yp) const {
	FT dummy(0);
	FT xs, ys, zs, zs_prev;
	m_p.projector().plane2Sphere(xp, yp, dummy, LIB_RATSS_NAMESPACE::SP_UPPER, xs, ys, zs);
	zs_prev = zs;
// 	std::cout << "zs=" << zs << std::endl;
// 	std::cout << "m_epsZ=" << m_epsZ << std::endl;
	
	uint32_t count = 0;
	while(zs < epsZ()) {
		xp /= 2;
		yp /= 2;
		zs_prev = zs;
		m_p.projector().plane2Sphere(xp, yp, dummy, LIB_RATSS_NAMESPACE::SP_UPPER, xs, ys, zs);
// 		std::cout << zs << std::endl;
		assert(zs_prev < zs);
		
		++count;
	}
// 	std::cout << "Generating an aux point took " << count << " rounds" << std::endl;
	return Point_3(xs, ys, zs);
}
//END private inital construction functions

#undef TMPL_CLS
#undef TMPL_HDR

} //end namespace LIB_DTS2_NAMESPACE



#endif
