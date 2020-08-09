#include <libdts2/constants.h>
#include <libdts2/Kernel_sp/Kernel_sp.h>
#include <libdts2/debug.h>

#include "TestBase.h"
#include "../vendor/libratss/common/generators.h"

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

void Kernel_sp_Test::setUp() {}

void Kernel_sp_Test::special() {
	dts2::detail::Point_sp_base pb(51406432, 16885681, 268435456, ratss::SP_DIM1_NEGATIVE);
	dts2::detail::Point_sp_base qb(102807927, 33763820, 536870912, ratss::SP_DIM1_NEGATIVE);
	dts2::detail::Point_sp_base rb(102783408, 33816407, 536870912, ratss::SP_DIM1_NEGATIVE);
	dts2::detail::Point_sp_base tb(102814795, 33777817, 536870912, ratss::SP_DIM1_NEGATIVE);
	{
		using Kernel = dts2::Kernel_sp;
		using Pred = Kernel::Side_of_oriented_circle_2;
		using Point = Kernel::Point_2;
		Kernel kernel;
		Pred pred = kernel.side_of_oriented_circle_2_object();
		Point p(pb), q(qb), r(rb), t(tb);
		auto result = pred(p, q, r, t);
		auto should = pred.MyBaseClass::operator()(p.point3(), q.point3(), r.point3(), t.point3());
		CPPUNIT_ASSERT_EQUAL(result, should);
	}
}

}} //end namespace ratss::tests
