#ifndef LIB_DTS2_POINT_SP_H
#define LIB_DTS2_POINT_SP_H
#pragma once

#include <libdts2/constants.h>
#include <libratss/CGAL/ExtendedInt64q.h>
#include <libratss/ProjectSN.h>

namespace LIB_DTS2_NAMESPACE {


template<typename T_BASE_TRAITS>
class Point_sp {
public:
	using MyBaseTrait = T_BASE_TRAITS;
	using base_type = int32_t;
	using BFT = MyBaseTrait::FT;
	using Point_3 = typename MyBaseTrait::Point_3;
public:
	Point_sp();
	Point_sp(const Point_sp & other);
	template<typename T_POINT_3>
	Point_sp(const T_POINT_3 & p);
	Point_sp(base_type num0, base_type num1, base_type pos, base_type exponent);
public:
	Point_sp & operator=(const Point_sp & other);
public:
	BFT x() const;
	BFT y() const;
	BFT z() const;
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
	int32_t m_nums[2];
	int32_t m_pos:4; //position on sphere, can be any of ratss::PositionOnSphere
	int32_t m_e:16; //denominator exponent denom = 1 << m_e
};

}//end namespace

//the definitions

namespace LIB_DTS2_NAMESPACE {
	
};

#endif