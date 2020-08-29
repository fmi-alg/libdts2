#include <libdts2/Kernel_sp/Point_sp.h>

#include <libratss/Conversion.h>
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
	if (Storage::dense_coding) {
		throw std::runtime_error("Unimplemented function");
	}
	else {
		if (sizeof(Types::numerator) != Storage::Size::NUMERATOR0) {
			bool sign = v<0;
			if (sign) {
				v = -v;
			}
			v <<= 0;
			if (sign) {
				v |= 1;
			}
		}
		::memmove(m_d.data()+Storage::Begin::NUMERATOR0, &v, Storage::Size::NUMERATOR0);
	}
}

void
Point_sp_base::set_numerator1(Types::numerator v) {
	if (Storage::dense_coding) {
		throw std::runtime_error("Unimplemented function");
	}
	else {
		if (sizeof(Types::numerator) != Storage::Size::NUMERATOR1) {
			bool sign = v<0;
			if (sign) {
				v = -v;
			}
			v <<= 0;
			if (sign) {
				v |= 1;
			}
		}
		::memmove(m_d.data()+Storage::Begin::NUMERATOR1, &v, Storage::Size::NUMERATOR1);
	}
}

namespace {
	std::size_t clz(uint64_t v) { return __builtin_clzl(v); }
	std::size_t clz(uint32_t v) { return __builtin_clz(v); }
	std::size_t clz(__uint128_t v) {
		uint64_t upper = v >> 64;
		uint64_t lower = v;
		if (upper) {
			return clz(upper);
		}
		else {
			return clz(lower)+64;
		}
	}
}

void
Point_sp_base::set_denominator(Types::denominator v) {
	if (v <= 0) {
		throw std::runtime_error("Point_sp_base::set_denominator: denominator has to be greater than 0");
	}
	if ((v & (v-1)) != 0) {
		throw std::runtime_error("Point_sp_base::set_denominator: denominator has to be a power of 2");
	}
	set_exponent((std::numeric_limits< std::make_unsigned<Types::numerator>::type >::digits-1)-clz(std::make_unsigned<Types::denominator>::type(v)));
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
		if (Storage::dense_coding) {
			throw std::runtime_error("Unimplemented function");
		}
		else {
			m_d.at(Storage::Begin::EXPONENT) = v;
		}
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
	m_d.at(Storage::Begin::POS) = v;
}

void
Point_sp_base::set_numerator0(mpz_class v) {
	set_numerator0( LIB_RATSS_NAMESPACE::convert<Types::numerator>(v) );
}

void
Point_sp_base::set_numerator1(mpz_class v) {
	set_numerator1( LIB_RATSS_NAMESPACE::convert<Types::numerator>(v) );
}

void
Point_sp_base::set_denominator(mpz_class v) {
	set_denominator( LIB_RATSS_NAMESPACE::convert<Types::denominator>(v) );
}

Point_sp_base::Types::numerator
Point_sp_base::numerator0() const {
	Types::numerator result = 0;
	if (Storage::dense_coding) {
		throw std::runtime_error("Unimplemented function");
	}
	else {
		::memmove(&result, m_d.data()+Storage::Begin::NUMERATOR0, Storage::Size::NUMERATOR0);
		if (sizeof(Types::numerator) != Storage::Size::NUMERATOR0) {
			bool sign = result & 0x1;
			result >>= 1;
			if (sign) {
				result = -result;
			}
		}
	}
	return result;
}

Point_sp_base::Types::numerator
Point_sp_base::numerator1() const {
	Types::numerator result = 0;
	if (Storage::dense_coding) {
		throw std::runtime_error("Unimplemented function");
	}
	else {
		::memmove(&result, m_d.data()+Storage::Begin::NUMERATOR1, Storage::Size::NUMERATOR1);
		if (sizeof(Types::numerator) != Storage::Size::NUMERATOR1) {
			bool sign = result & 0x1;
			result >>= 1;
			if (sign) {
				result = -result;
			}
		}
	}
	return result;
}

Point_sp_base::Types::denominator
Point_sp_base::denominator() const {
	assert(exponent() < Digits::DENOMINATOR);
	return static_cast<Types::denominator>(1) << exponent();
}

LIB_RATSS_NAMESPACE::PositionOnSphere
Point_sp_base::pos() const {
	if (Storage::dense_coding) {
		throw std::runtime_error("Unimplemented function");
	}
	else {
		return static_cast<LIB_RATSS_NAMESPACE::PositionOnSphere>(m_d.at(Storage::Begin::POS));
	}
}

uint8_t
Point_sp_base::exponent() const {
	if (fixed_exponent) {
		return fixed_exponent;
	}
	else {
		if (Storage::dense_coding) {
			throw std::runtime_error("Unimplemented function");
		}
		else {
			return m_d.at(Storage::Begin::EXPONENT);
		}
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
