#pragma once
#ifndef LIB_DTS2_POINT_SP_H
#define LIB_DTS2_POINT_SP_H

#include <libdts2/constants.h>
#include <libratss/CGAL/ExtendedInt64q.h>
#include <libratss/ProjectS2.h>

namespace LIB_DTS2_NAMESPACE {


template<typename T_BASE_TRAITS>
class Point_sp {
public:
	using MyBaseTrait = T_BASE_TRAITS;
	using base_type = int32_t;
	using BFT = MyBaseTrait::FT;
	using Point_3 = typename MyBaseTrait::Point_3;
	static constexpr int snap_bits = 30;
public:
	Point_sp();
	Point_sp(const Point_sp & other);
	template<typename T_POINT_3>
	Point_sp(const T_POINT_3 & p);
	Point_sp(base_type num0, base_type num1, base_type denominator, base_type pos);
public:
	Point_sp & operator=(const Point_sp & other);
public:
	BFT x() const;
	BFT y() const;
	BFT z() const;
	Point_3 point3() const;
	operator Point_3() const;
public:
	base_type numerator0() const;
	base_type numerator1() const;
	///if this is an auxiliary point, then this will overflow
	base_type denominator() const;
	///@return position on Sphere
	int32_t pos() const;
	int32_t exponent() const;
	bool is_auxiliary() const;
public:
	static Point_sp make_auxiliary();
private:
	struct Nums {
		int32_t nums[2];
		int32_t dummy;
		Nums(base_type num0, base_type num1) : nums({num0, num1}), dummy(0) {}
	};
	struct Meta {
		int32_t dummy[2];
		int32_t pos:4; //position on sphere, can be any of ratss::PositionOnSphere
		int32_t e:16;  //denominator exponent denom = 1 << m_e
		int32_t reserved:12;
	};
	struct Raw {
		int32_t d[3];
		Raw() : d({0,0,0}) {}
	};
	union Data {
		Nums nums;
		Meta meta;
		Raw raw;
		Data() : raw() {}
		Data(base_type num0, base_type num1) : nums(num0, num1) {}
	};
private:
	Data m_d;
};

}//end namespace

//the definitions

namespace LIB_DTS2_NAMESPACE {

#define PTSP_TMP_PRMS template<typename T_BASE_TRAITS>
#define PTSP_CLS_NAME Point_sp<T_BASE_TRAITS>

PTSP_TMP_PRMS
PTSP_CLS_NAME::Point_sp() :
m_d()
{}


PTSP_TMP_PRMS
PTSP_CLS_NAME::Point_sp(const Point_sp & other) :
m_d(other.m_d)
{}


PTSP_TMP_PRMS
template<typename T_POINT_3>
PTSP_CLS_NAME::Point_sp(const T_POINT_3 & p) {
	using FT = T_POINT_3::FT;
	
	
	LIB_RATSS_NAMESPACE::ProjectS2 proj;
	proj.sphere2Plane();
}

PTSP_TMP_PRMS
PTSP_CLS_NAME::Point_sp(base_type num0, base_type num1, base_type denominator, base_type pos) :
m_d(num0, num1)
{
	using std::abs;
	//based on https://stackoverflow.com/a/108360
	assert(denominator > 0 && unsigned_base_type(denominator) && (denominator & (denominator-1)) == 0);
	assert(abs(num0) <= denominator && abs(num1) <= denominator);
	m_d.meta.pos = pos;
	m_d.meta.e = 31-__builtin_clz(denominator);
}

#undef PTSP_TMP_PRMS
#undef PTSP_CLS_NAME

};

#endif
