#pragma once

#ifndef LIB_DTS2_POINT_SP_H
#define LIB_DTS2_POINT_SP_H

#include <libdts2/constants.h>
#include <libratss/ProjectS2.h>

#include <array>
#include <stdint.h>
#include <string.h>
#include <assert.h>

namespace LIB_DTS2_NAMESPACE {
namespace detail {

class Point_sp_base {
public:
	using base_type = int64_t;
	using unsigned_base_type = uint64_t;
	//TODO: extend this to 42 bits for numerators
	struct BitSizes {
		static constexpr int NUMERATOR0=32;
		static constexpr int NUMERATOR1=32;
		static constexpr int DENOMINATOR=7;
		static constexpr int POS=3;
	};
public:
	Point_sp_base();
	Point_sp_base(Point_sp_base const &) = default;
	Point_sp_base(base_type _num0, base_type _num1, unsigned_base_type _den, base_type _pos);
	virtual ~Point_sp_base();
	Point_sp_base & operator=(Point_sp_base const &) = default;
public:
	void set_numerator0(base_type v);
	void set_numerator1(base_type v);
	void set_denominator(unsigned_base_type v);
	void set_exponent(uint8_t v);
	void set_pos(int v);
	void set_numerator0(mpz_class v);
	void set_numerator1(mpz_class v);
	void set_denominator(mpz_class v);
public:
	base_type numerator0() const;
	base_type numerator1() const;
	unsigned_base_type denominator() const;
	///@return position on Sphere
	LIB_RATSS_NAMESPACE::PositionOnSphere pos() const;
	uint8_t exponent() const;
public:
	bool operator!=(Point_sp_base const & other) const;
	bool operator==(Point_sp_base const & other) const;
private:
	//Coordinates in the plane are of the form
	//x = num/2^s with s < 128
	//Thus we use 7 Bits to encode s
	//Additionally we need 3 Bits to encode the 6 position on the sphere
	//The remaining 84 Bits are used for the coordinates
	std::array<char, 12> m_d;
};
	
} //end namespace detail

template<typename T_BASE_TRAITS>
class Point_sp: public detail::Point_sp_base {
public:
	using MyParent = detail::Point_sp_base;
	using MyBaseTrait = T_BASE_TRAITS;
	using FT = typename MyBaseTrait::FT;
	using Point_3 = typename MyBaseTrait::Point_3;
public:
	Point_sp();
	Point_sp(MyParent const & v);
	Point_sp(base_type _num0, base_type _num1, unsigned_base_type _den, base_type _pos);
	Point_sp(FT const & x, FT const & y, FT const & z);
	Point_sp(Point_3 const & v);
	Point_sp(Point_sp const & other);
public:
	Point_sp & operator=(Point_sp const & other) = default;
public:
	Point_3 point3() const;
public:
	FT x() const { return point3().x(); }
	FT y() const { return point3().y(); }
	FT z() const { return point3().z(); }
public:
	operator Point_3() const { return point3(); }
};

template<typename K1, typename K2>
bool operator==(Point_sp<K1> const & a, Point_sp<K2> const & b) {
	return static_cast<detail::Point_sp_base const &>(a) == static_cast<detail::Point_sp_base const &>(b);
}

template<typename K1, typename K2>
bool operator!=(Point_sp<K1> const & a, Point_sp<K2> const & b) {
	return static_cast<detail::Point_sp_base const &>(a) != static_cast<detail::Point_sp_base const &>(b);
}

}//end namespace

//the definitions

namespace LIB_DTS2_NAMESPACE {

#define PTSP_TMP_PRMS template<typename T_BASE_TRAITS>
#define PTSP_CLS_NAME Point_sp<T_BASE_TRAITS>

PTSP_TMP_PRMS
PTSP_CLS_NAME::Point_sp() {}

PTSP_TMP_PRMS
PTSP_CLS_NAME::Point_sp(const Point_sp & other) :
MyParent(other)
{}

PTSP_TMP_PRMS
PTSP_CLS_NAME::Point_sp(Point_3 const & p) :
Point_sp(p.x(), p.y(), p.z())
{
	assert(point3() == p);
}

PTSP_TMP_PRMS
PTSP_CLS_NAME::Point_sp(FT const & x, FT const & y, FT const & z) {
	if (x == LIB_DTS2_ORIGIN_X && y == LIB_DTS2_ORIGIN_Y && z == LIB_DTS2_ORIGIN_Z) { //Hack for LIB_DTS2_ORIGIN 
		set_pos(0);
		return;
	}
	assert(x*x + y*y + z*z == FT(1));
	LIB_RATSS_NAMESPACE::ProjectS2 proj;
	using LIB_RATSS_NAMESPACE::convert;
	std::array<mpq_class, 3> pp;
	int _pos = proj.sphere2Plane(
								 convert<mpq_class>(x),
								 convert<mpq_class>(y),
								 convert<mpq_class>(z),
								 pp[0], pp[1], pp[2]
								);
	this->set_pos(_pos);
	//get maximumum of all denominators
	using std::max;
	mpz_class mden = max(pp[0].get_den(), max(pp[1].get_den(), pp[2].get_den()));
	
	set_denominator(mden);
	
	switch(abs(_pos)) {
	case 1:
		set_numerator0(pp[1].get_num()*(mden/pp[1].get_den()));
		set_numerator1(pp[2].get_num()*(mden/pp[2].get_den()));
		break;
	case 2:
		set_numerator0(pp[0].get_num()*(mden/pp[0].get_den()));
		set_numerator1(pp[2].get_num()*(mden/pp[2].get_den()));
		break;
	case 3:
		set_numerator0(pp[0].get_num()*(mden/pp[0].get_den()));
		set_numerator1(pp[1].get_num()*(mden/pp[1].get_den()));
		break;
	default:
		throw std::runtime_error("Point_sp: Position is out of bounds");
	}
	assert(Point_3(x,y,z) == point3());
}

PTSP_TMP_PRMS
PTSP_CLS_NAME::Point_sp(base_type _num0, base_type _num1, unsigned_base_type _den, base_type _pos) :
MyParent(_num0, _num1, _den, _pos)
{}

PTSP_TMP_PRMS
PTSP_CLS_NAME::Point_sp(MyParent const & v) :
MyParent(v)
{}

PTSP_TMP_PRMS
typename PTSP_CLS_NAME::Point_3
PTSP_CLS_NAME::point3() const {
	if (pos() == 0) {  //Hack for LIB_DTS2_ORIGIN 
		return LIB_DTS2_ORIGIN;
	}
	
	LIB_RATSS_NAMESPACE::ProjectS2 proj;
	mpq_class xp, yp, zp;
	mpq_class xs, ys, zs;
	mpq_class den = denominator();
	
	switch (std::abs(this->pos())) {
	case 1:
		xp = 0;
		yp = mpq_class(numerator0())/den;
		zp = mpq_class(numerator1())/den;
		break;
	case 2:
		xp = mpq_class(numerator0())/den;
		yp = 0;
		zp = mpq_class(numerator1())/den;
		break;
	case 3:
		xp = mpq_class(numerator0())/den;
		yp = mpq_class(numerator1())/den;
		zp = 0;
		break;
	default:
		throw std::runtime_error("Point_sp: Position is out of bounds");
	};
	proj.plane2Sphere(xp, yp, zp, pos(), xs, ys, zs);
	
	using LIB_RATSS_NAMESPACE::convert;
	return Point_3( convert<FT>(xs), convert<FT>(ys), convert<FT>(zs) );
}

#undef PTSP_TMP_PRMS
#undef PTSP_CLS_NAME

};

#endif
