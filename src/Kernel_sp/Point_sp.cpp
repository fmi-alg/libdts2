#include <libdts2/Kernel_sp/Point_sp.h>

#include <limits>

namespace LIB_DTS2_NAMESPACE {
namespace detail {

Point_sp_base::Point_sp_base() {
	m_d.fill(0);
}

Point_sp_base::Point_sp_base(Types::numerator _num0, Types::numerator _num1, Types::denominator _den, Types::pos _pos) {
	set_numerator0(_num0);
	set_numerator1(_num1);
	set_denominator(_den);
	set_pos(_pos);
}

Point_sp_base::~Point_sp_base() {}

void
Point_sp_base::set_numerator0(Types::numerator v) {
	assert(int32_t(v) == v);
	int32_t tv(v);
	::memmove(m_d.data(), &tv, sizeof(tv));
}

void
Point_sp_base::set_numerator1(Types::numerator v) {
	assert(int32_t(v) == v);
	int32_t tv(v);
	::memmove(m_d.data()+4, &tv, sizeof(tv));
}

namespace {
	std::size_t clz(uint64_t v) { return __builtin_clzl(v); }
	std::size_t clz(uint32_t v) { return __builtin_clz(v); }
}

void
Point_sp_base::set_denominator(Types::denominator v) {
	if (v <= 0) {
		throw std::runtime_error("Point_sp_base::set_denominator: denominator has to be greater than 0");
	}
	if ((v & (v-1)) != 0) {
		throw std::runtime_error("Point_sp_base::set_denominator: denominator has to be a power of 2");
	}
	set_exponent((std::numeric_limits<Types::numerator>::digits-1)-clz(std::make_unsigned<Types::denominator>::type(v)));
	assert(denominator() == v);
}

void 
Point_sp_base::set_exponent(Types::exponent v) {
	if (fixed_exponent) {
		if (v != fixed_exponent) {
			throw std::runtime_error("Point_sp_base: fixed exponent does not match set exponent");
		}
	}
	else {
		m_d.at(9) = v;
	}
}

//map is
//0 -> 0
//+1 -> 1
//-1 -> 2
//+2 -> 3
//-2 -> 4
//+3 -> 5
//-3 -> 6
//This boils down to
//v==0 -> v
//v & 0x1 -> (v+1)/2
//!(v & 0x1) -> -v/2

void 
Point_sp_base::set_pos(Types::pos v) {
	assert(Types::pos(v) == v);
	m_d.at(8) = v;
}

void
Point_sp_base::set_numerator0(mpz_class v) {
	assert(v.fits_slong_p());
	assert(Types::numerator(v.get_si()) == v.get_si());
	set_numerator0(v.get_si());
}

void
Point_sp_base::set_numerator1(mpz_class v) {
	assert(v.fits_slong_p());
	assert(Types::numerator(v.get_si()) == v.get_si());
	set_numerator1(v.get_si());
}

void
Point_sp_base::set_denominator(mpz_class v) {
	assert(v.fits_slong_p());
	assert(Types::denominator(v.get_si()) == v.fits_slong_p());
	set_denominator(Types::denominator(v.get_si()));
}

Point_sp_base::Types::numerator
Point_sp_base::numerator0() const {
	static_assert(std::is_integral<Types::numerator>::value);
	Types::numerator tmp = 0;
	::memmove(&tmp, m_d.data(), sizeof(tmp));
	return tmp;
}

Point_sp_base::Types::numerator
Point_sp_base::numerator1() const {
	static_assert(std::is_integral<Types::numerator>::value);
	Types::numerator tmp = 0;
	::memmove(&tmp, m_d.data()+sizeof(Types::numerator), sizeof(tmp));
	return tmp;
}

Point_sp_base::Types::denominator
Point_sp_base::denominator() const {
	assert(exponent() < Digits::DENOMINATOR);
	return static_cast<Types::denominator>(1) << exponent();
}

LIB_RATSS_NAMESPACE::PositionOnSphere
Point_sp_base::pos() const {
	return static_cast<LIB_RATSS_NAMESPACE::PositionOnSphere>(m_d.at(8));
}

uint8_t
Point_sp_base::exponent() const {
	if (fixed_exponent) {
		return fixed_exponent;
	}
	else {
		return m_d.at(9);
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
