#include <libdts2/constants.h>
#include <libdts2/Kernel_sp/Kernel_sp.h>
#include <libdts2/debug.h>

#include "TestBase.h"
#include "../vendor/libratss/common/generators.h"

#include <type_traits>

namespace LIB_DTS2_NAMESPACE {
namespace tests {

class Kernel_sp_Test: public TestBase {
CPPUNIT_TEST_SUITE( Kernel_sp_Test );
CPPUNIT_TEST( special );
CPPUNIT_TEST_SUITE_END();
public:
	static std::size_t num_random_test_points;
public:
	Kernel_sp_Test() {}
public:
	virtual void setUp() override;
public:
	void special();
};

std::size_t Kernel_sp_Test::num_random_test_points;


}} // end namespace dts2::tests
int main(int argc, char ** argv) {
	dts2::init_interactive_debuging();
	dts2::tests::TestBase::init(argc, argv);
	dts2::tests::Kernel_sp_Test::num_random_test_points = 10000;
	srand( 0 );
	CppUnit::TextUi::TestRunner runner;
	runner.addTest(  LIB_DTS2_NAMESPACE::tests::Kernel_sp_Test::suite() );
	bool ok = runner.run();
	return ok ? 0 : 1;
}

namespace LIB_DTS2_NAMESPACE {
namespace tests {
	
struct CartesianRational {
	std::array<mpq_class, 3> coords;
	CartesianRational() {}
	CartesianRational(mpq_class const & x, mpq_class const & y, mpq_class const & z) :
	coords({x, y, z})
	{}
	mpq_class & operator[](std::size_t i) { return coords[i]; }
	mpq_class const & operator[](std::size_t i) const { return coords[i]; }
	CartesianRational snap() const {
		ratss::ProjectSN proj;
		CartesianRational pp, out;
		auto pos = proj.sphere2Plane(coords.begin(), coords.end(), pp.coords.begin());
		for(auto & x : pp.coords) {
			x = proj.calc().snap(x, ratss::ST_FX, 30);
		}
		proj.plane2Sphere(pp.coords.begin(), pp.coords.end(), pos, out.coords.begin());
		return out;
	}
	
	template<typename T_POINT, typename FT = typename T_POINT::FT>
	T_POINT point() const {
		return T_POINT(
			ratss::convert<FT>(coords[0]),
			ratss::convert<FT>(coords[1]),
			ratss::convert<FT>(coords[2])
		);
	}
	ratss::PositionOnSphere pos() const {
		ratss::ProjectSN proj;
		return proj.positionOnSphere(coords.begin(), coords.end());
	}
	
	void cycle_coords() {
		using std::swap;
		swap(coords[0], coords[2]);
		swap(coords[1], coords[2]);
	}
};

std::ostream & operator<<(std::ostream & out, CartesianRational const & v) {
	out << "(" << v.coords[0] << ", ";
	out << v.coords[1] << ", ";
	out << v.coords[2] << ")";
	return out;
}

void Kernel_sp_Test::setUp() {}

void Kernel_sp_Test::special() {
	using CGK = CGAL::Exact_predicates_exact_constructions_kernel;
	using PointCGK = CGK::Point_3;
	using Kernel = dts2::Kernel_sp;
	using Pred = Kernel::Side_of_oriented_circle_2;
	using Point = Kernel::Point_2;

	CGK cgk;
	auto pred_cgk = cgk.orientation_3_object();
	
	Kernel kernel;
	Pred pred = kernel.side_of_oriented_circle_2_object();
	
	CartesianRational pc(1, -0.25, 0);
	CartesianRational qc(1, 0.25, 0);
	CartesianRational rc(1, 0, 0.25);
	CartesianRational tic(1, 0, 0);
	CartesianRational toc(1, 0, 0.25+0.125);

	for(int i(1); i <= 3; ++i) {
		CartesianRational ps = pc.snap();
		CartesianRational qs = qc.snap();
		CartesianRational rs = rc.snap();
		CartesianRational tis = tic.snap();
		CartesianRational tos = toc.snap();
		
		std::cout << std::endl;
	#define PP(__NAME) std::cout << #__NAME "c" << __NAME ## c << " -> " #__NAME "s" << __NAME ## s << " -> " << (__NAME ## s).point<Point>() << std::endl;
		PP(p);
		PP(q);
		PP(r);
		PP(ti);
		PP(to);
	#undef PP
		
		CPPUNIT_ASSERT_EQUAL_MESSAGE("ps.pos", i, (int)ps.pos());
		CPPUNIT_ASSERT_EQUAL_MESSAGE("qs.pos",i, (int)qs.pos());
		CPPUNIT_ASSERT_EQUAL_MESSAGE("rs.pos",i, (int)rs.pos());
		CPPUNIT_ASSERT_EQUAL_MESSAGE("tis.pos",i, (int)tis.pos());
		CPPUNIT_ASSERT_EQUAL_MESSAGE("tos.pos",i, (int)tos.pos());
		
		PointCGK p_cgk = ps.point<PointCGK, CGK::FT>();
		PointCGK q_cgk = qs.point<PointCGK, CGK::FT>();
		PointCGK r_cgk = rs.point<PointCGK, CGK::FT>();
		PointCGK ti_cgk = tis.point<PointCGK, CGK::FT>();
		PointCGK to_cgk = tos.point<PointCGK, CGK::FT>();
		
		Point p = ps.point<Point>();
		Point q = qs.point<Point>();
		Point r = rs.point<Point>();
		Point ti = tis.point<Point>();
		Point to = tos.point<Point>();

		auto should_i_cgk = pred_cgk(p_cgk, q_cgk, r_cgk, ti_cgk);
		auto should_o_cgk = pred_cgk(p_cgk, q_cgk, r_cgk, to_cgk);
		
		auto result_i = pred(p, q, r, ti);
		auto should_i = pred.MyBaseClass::operator()(p.point3(), q.point3(), r.point3(), ti.point3());
		
		auto result_o = pred(p, q, r, to);
		auto should_o = pred.MyBaseClass::operator()(p.point3(), q.point3(), r.point3(), to.point3());
		
		CPPUNIT_ASSERT_EQUAL_MESSAGE("should_i_cgk", CGAL::ON_POSITIVE_SIDE, should_i_cgk);
		CPPUNIT_ASSERT_EQUAL_MESSAGE("should_o_cgk", CGAL::ON_NEGATIVE_SIDE, should_o_cgk);
		
		CPPUNIT_ASSERT_EQUAL_MESSAGE("should_i", CGAL::ON_POSITIVE_SIDE, should_i);
		CPPUNIT_ASSERT_EQUAL_MESSAGE("should_o", CGAL::ON_NEGATIVE_SIDE, should_o);
		
		CPPUNIT_ASSERT_EQUAL_MESSAGE("should_i_cgk == should_i", should_i_cgk, should_i);
		CPPUNIT_ASSERT_EQUAL_MESSAGE("should_o_cgk == should_o", should_o_cgk, should_o);
		
		CPPUNIT_ASSERT_EQUAL_MESSAGE("result_i", result_i, should_i);
		CPPUNIT_ASSERT_EQUAL_MESSAGE("result_o", result_o, should_o);
		
		pc.cycle_coords();
		qc.cycle_coords();
		rc.cycle_coords();
		tic.cycle_coords();
		toc.cycle_coords();
	}
}

}} //end namespace ratss::tests
