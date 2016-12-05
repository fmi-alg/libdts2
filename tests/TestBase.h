#ifndef DTS2_TESTS_TEST_BASE_H
#define DTS2_TESTS_TEST_BASE_H
#pragma once
#include <cppunit/ui/text/TestRunner.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/Asserter.h>
#include <cppunit/TestAssert.h>
#include <cppunit/TestResult.h>

#include <libdts2/constants.h>

namespace LIB_DTS2_NAMESPACE {
namespace tests {

class TestBase: public CppUnit::TestFixture {
// CPPUNIT_TEST_SUITE( Test );
// CPPUNIT_TEST( test );
// CPPUNIT_TEST_SUITE_END();
public:
	TestBase();
	virtual ~TestBase();
public:
	static void init(int argc, char ** argv);
private:
	static int argc;
	static char ** argv;
	
};

}} //end namespace dts2::tests

/*
int main(int argc, char ** argv) {
	dts2::tests::TestBase::init(argc, argv);
	srand( 0 );
	CppUnit::TextUi::TestRunner runner;
	runner.addTest(  Test::suite() );
	bool ok = runner.run();
	return ok ? 0 : 1;
}
*/

#endif