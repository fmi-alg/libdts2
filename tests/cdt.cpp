#include <libdts2/constants.h>
#include <libdts2/Constrained_delaunay_triangulation_s2.h>

#include "TestBase.h"
#include <vendor/libratss/common/generators.h>

namespace LIB_DTS2_NAMESPACE {
namespace tests {

class CDTTest: public TestBase {
CPPUNIT_TEST_SUITE( CDTTest );
CPPUNIT_TEST( valid );
CPPUNIT_TEST_SUITE_END();
public:
	static std::size_t num_random_test_points;
private:
	Constrained_Delaunay_triangulation_no_intersections_with_info_s2<void, void> m_cdt;
public:
	CDTTest() : m_cdt(53) {}
public:
	virtual void setUp() override;
public:
	void valid();
};

std::size_t CDTTest::num_random_test_points;


}} // end namespace dts2::tests
int main(int argc, char ** argv) {
	dts2::tests::TestBase::init(argc, argv);
	dts2::tests::CDTTest::num_random_test_points = 10000;
	srand( 0 );
	CppUnit::TextUi::TestRunner runner;
	runner.addTest(  LIB_DTS2_NAMESPACE::tests::CDTTest::suite() );
	bool ok = runner.run();
	return ok ? 0 : 1;
}

namespace LIB_DTS2_NAMESPACE {
namespace tests {

void CDTTest::setUp() {
	std::vector<ratss::SphericalCoord> coords = ratss::getRandomPolarPoints(CDTTest::num_random_test_points);
	m_cdt.insert(coords.begin(), coords.end());
}

void CDTTest::valid() {
	CPPUNIT_ASSERT(m_cdt.is_valid());
}

}} //end namespace ratss::tests