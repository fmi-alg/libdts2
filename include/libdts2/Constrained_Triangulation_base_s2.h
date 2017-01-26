#ifndef LIB_DTS2_CONSTRAINED_TRIANGULATION_BASE_S2_H
#define LIB_DTS2_CONSTRAINED_TRIANGULATION_BASE_S2_H
#pragma once
#include <libdts2/constants.h>
#include <libdts2/Triangulation_base_s2.h>

namespace LIB_DTS2_NAMESPACE {

template<typename T_BASE_TRIANGULATION>
class Constrained_triangulation_base_s2: public Triangulation_base_s2<T_BASE_TRIANGULATION> {
public:
	using MyBaseClass = Triangulation_base_s2<T_BASE_TRIANGULATION>;
public:
	using Geom_traits = typename MyBaseClass::Geom_traits;
	using Intersection_tag = typename T_BASE_TRIANGULATION::Intersection_tag;
	using Itag = Intersection_tag;
public:
	Constrained_triangulation_base_s2(int significands);
	Constrained_triangulation_base_s2(const Geom_traits & traits);
	Constrained_triangulation_base_s2(Constrained_triangulation_base_s2 && other);
	Constrained_triangulation_base_s2(const Constrained_triangulation_base_s2 & other) = delete;
	Constrained_triangulation_base_s2 & operator=(Constrained_triangulation_base_s2 && other);
};

}//end namespace LIB_DTS2_NAMESPACE

//definition
namespace LIB_DTS2_NAMESPACE {
#define TMPL_HDR template<typename T_BASE_TRIANGULATION>
#define TMPL_CLS Constrained_triangulation_base_s2<T_BASE_TRIANGULATION>

//BEGIN constructors
TMPL_HDR
TMPL_CLS::Constrained_triangulation_base_s2(int significands) :
MyBaseClass(significands)
{}


TMPL_HDR
TMPL_CLS::Constrained_triangulation_base_s2(const Geom_traits & traits) :
MyBaseClass(traits)
{}

TMPL_HDR
TMPL_CLS::Constrained_triangulation_base_s2(Constrained_triangulation_base_s2 && other) :
MyBaseClass(std::move(other))
{}

TMPL_HDR
TMPL_CLS & 
TMPL_CLS::operator=(TMPL_CLS && other) {
	MyBaseClass::operator=(std::move(other));
	return *this;
}
//END constructors
#undef TMPL_CLS
#undef TMPL_HDR

}//end namespace LIB_DTS2_NAMESPACE

#endif