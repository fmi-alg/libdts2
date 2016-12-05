#include <libdts2/constants.h>
#include <libdts2/Delaunay_triangulation_s2.h>

#include "TestBase.h"
#include <vendor/libratss/common/generators.h>

namespace LIB_DTS2_NAMESPACE {
namespace tests {

class DTTest: public TestBase {
CPPUNIT_TEST_SUITE( DTTest );
CPPUNIT_TEST( valid );
CPPUNIT_TEST_SUITE_END();
public:
	static std::size_t num_random_test_points;
private:
	Delaunay_triangulation_with_info_s2<void, void> m_cdt;
public:
	DTTest() : m_cdt(53) {}
public:
	virtual void setUp() override;
public:
	void valid();
};

std::size_t DTTest::num_random_test_points;


}} // end namespace dts2::tests
int main(int argc, char ** argv) {
	dts2::tests::TestBase::init(argc, argv);
	dts2::tests::DTTest::num_random_test_points = 10000;
	srand( 0 );
	CppUnit::TextUi::TestRunner runner;
	runner.addTest(  LIB_DTS2_NAMESPACE::tests::DTTest::suite() );
	bool ok = runner.run();
	return ok ? 0 : 1;
}

namespace LIB_DTS2_NAMESPACE {
namespace tests {

void DTTest::setUp() {
	std::vector<ratss::SphericalCoord> coords = ratss::getRandomPolarPoints(DTTest::num_random_test_points);
	m_cdt.insert(coords.begin(), coords.end());
}

void DTTest::valid() {
	CPPUNIT_ASSERT(m_cdt.is_valid());
}

}} //end namespace ratss::tests