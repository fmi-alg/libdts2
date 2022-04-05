#include <libdts2/debug.h>
#include <libratss/debug.h>
#include <libratss/Conversion.h>

extern "C" {
	
void dts2_debug_print_point_sp(LIB_DTS2_NAMESPACE::Constrained_delaunay_triangulation_with_inexact_intersections_traits_s2_sp::Point_2 const & p) {
	using ratss::convert;
	std::cerr << "Point_sp(num0=" << p.numerator0() << ", num1=" << p.numerator1() << ", den=" << p.denominator() << ", pos=" << p.pos() << ")=";
	std::cerr << "Point_3(" << std::endl;
	std::cerr << "x=";
	ratss_debug_print_mpq_class( convert<mpq_class>(p.x()) );
	std::cerr << "y=";
	ratss_debug_print_mpq_class( convert<mpq_class>(p.y()) );
	std::cerr << "z=";
	ratss_debug_print_mpq_class( convert<mpq_class>(p.z()) );
	std::cerr << ")" << std::endl;
}


void dts2_debug_print_kernel_sp_point(LIB_DTS2_NAMESPACE::Kernel_sp::Point_2 const & p) {
	using ratss::convert;
	std::cerr << "Point_sp(num0=" << p.numerator0() << ", num1=" << p.numerator1() << ", den=" << p.denominator() << ", pos=" << p.pos() << ")=";
	std::cerr << "Point_3(" << std::endl;
	std::cerr << "x=";
	ratss_debug_print_mpq_class( convert<mpq_class>(p.x()) );
	std::cerr << "y=";
	ratss_debug_print_mpq_class( convert<mpq_class>(p.y()) );
	std::cerr << "z=";
	ratss_debug_print_mpq_class( convert<mpq_class>(p.z()) );
	std::cerr << ")" << std::endl;
}

void dts2_debug_print_kernel_sp64_point(LIB_DTS2_NAMESPACE::Kernel_sp_64::Point_2 const & p) {
	using ratss::convert;
	std::cerr << "Point_sp(num0=" << p.numerator0() << ", num1=" << p.numerator1() << ", den=" << p.denominator() << ", pos=" << p.pos() << ")=";
	std::cerr << "Point_3(" << std::endl;
	std::cerr << "x=";
	ratss_debug_print_mpq_class( convert<mpq_class>(p.x()) );
	std::cerr << "y=";
	ratss_debug_print_mpq_class( convert<mpq_class>(p.y()) );
	std::cerr << "z=";
	ratss_debug_print_mpq_class( convert<mpq_class>(p.z()) );
	std::cerr << ")" << std::endl;
}

}

namespace LIB_DTS2_NAMESPACE {
	
#define PRINT_FN_ADDR(__F) std::cout << #__F ": " << std::size_t(&__F) << std::endl;
void init_interactive_debuging() {
	ratss::init_interactive_debuging();
	std::cout << "Initializing interactive debugging functionality of dts2:\n";
	PRINT_FN_ADDR(dts2_debug_print_point_sp)
	PRINT_FN_ADDR(dts2_debug_print_kernel_sp_point)
	PRINT_FN_ADDR(dts2_debug_print_kernel_sp64_point)
	std::cout << std::endl;
}
#undef PRINT_FN_ADDR

}
