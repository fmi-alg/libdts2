#include "TestBase.h"

namespace dts2 {
namespace tests {

int TestBase::argc = 0;
char ** TestBase::argv = 0;

TestBase::TestBase() {
// 	if (!argc) {
// 		throw std::runtime_error("sserialize::tests::TestBase:: need to call init first!");
// 	}
}

TestBase::~TestBase() {}

void TestBase::init(int argc, char ** argv) {
	TestBase::argc = argc;
	TestBase::argv = argv;
}


}} //end namespace dts2:.tests