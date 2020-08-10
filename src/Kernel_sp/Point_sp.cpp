#include <libdts2/Kernel_sp/Point_sp.h>

#include <limits>

namespace LIB_DTS2_NAMESPACE {
namespace detail {

Point_sp_base::Point_sp_base() {
	m_d.fill(0);
}

Point_sp_base::Point_sp_base(base_type _num0, base_type _num1, base_type _den, base_type _pos) {
	set_numerator0(_num0);
	set_numerator1(_num1);
	set_denominator(_den);
	set_pos(_pos);
}

Point_sp_base::~Point_sp_base() {}

void
Point_sp_base::set_numerator0(base_type v) {
	assert(int32_t(v) == v);
	int32_t tv(v);
	::memmove(m_d.data(), &tv, sizeof(tv));
}

void
Point_sp_base::set_numerator1(base_type v) {
	assert(int32_t(v) == v);
	int32_t tv(v);
	::memmove(m_d.data()+4, &tv, sizeof(tv));
}

void
Point_sp_base::set_denominator(base_type v) {
	if (v <= 0) {
		throw std::runtime_error("Point_sp_base::set_denominator: denominator has to be greater than 0");
	}
	if ((v & (v-1)) != 0) {
		throw std::runtime_error("Point_sp_base::set_denominator: denominator has to be a power of 2");
	}
	set_exponent((std::numeric_limits<unsigned_base_type>::digits-1)-__builtin_clzl(unsigned_base_type(v)));
	assert(denominator() == v);
}

void 
Point_sp_base::set_exponent(uint8_t v) {
	if (fixed_exponent) {
		if (v != fixed_exponent) {
			throw std::runtime_error("Point_sp_base: fixed exponent does not match set exponent");
		}
	}
	else {
		m_d[9] = v;
	}
}

void 
Point_sp_base::set_pos(int v) {
	assert(char(v) == v);
	m_d[8] = v;
}

void
Point_sp_base::set_numerator0(mpz_class v) {
	assert(v.fits_slong_p());
	set_numerator0(v.get_si());
}

void
Point_sp_base::set_numerator1(mpz_class v) {
	assert(v.fits_slong_p());
	set_numerator1(v.get_si());
}

void
Point_sp_base::set_denominator(mpz_class v) {
	assert(v.fits_slong_p());
	set_denominator(v.get_si());
}

Point_sp_base::base_type
Point_sp_base::numerator0() const {
	int32_t tmp = 0;
	::memmove(&tmp, m_d.data(), sizeof(tmp));
	return tmp;
}

Point_sp_base::base_type
Point_sp_base::numerator1() const {
	int32_t tmp = 0;
	::memmove(&tmp, m_d.data()+4, sizeof(tmp));
	return tmp;
}

Point_sp_base::base_type
Point_sp_base::denominator() const {
	return static_cast<unsigned_base_type>(1) << exponent();
}

LIB_RATSS_NAMESPACE::PositionOnSphere
Point_sp_base::pos() const {
	return static_cast<LIB_RATSS_NAMESPACE::PositionOnSphere>(m_d[8]);
}

uint8_t
Point_sp_base::exponent() const {
	if (fixed_exponent) {
		return fixed_exponent;
	}
	else {
		return m_d[9];
	}
}


bool
Point_sp_base::operator!=(Point_sp_base const & other) const {
	return m_d != other.m_d;
}

bool
Point_sp_base::operator==(Point_sp_base const & other) const {
	return m_d == other.m_d;
}

std::ostream & operator<<(std::ostream & out, Point_sp_base const & v) {
	out << "Point_sp"
		<< "(num0=" << v.numerator0()
		<< ", num1=" << v.numerator1()
		<< ", den=" << v.denominator()
		<< ", e=" << int(v.exponent())
		<< ", pos=" << int(v.pos())
		<< ")";
	return out;
}

}//end namesspace detail
}//end namespace LIB_DTS2_NAMESPACE
