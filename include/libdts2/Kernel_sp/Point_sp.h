#pragma once

#ifndef LIB_DTS2_POINT_SP_H
#define LIB_DTS2_POINT_SP_H

#include <libdts2/constants.h>
#include <libdts2/Kernel_sp/internal.h>
#include <libratss/ProjectS2.h>

#include <array>
#include <stdint.h>
#include <string.h>
#include <assert.h>

namespace LIB_DTS2_NAMESPACE {
namespace detail {
	
//TODO:Check endianness
	
struct Point_sp_cfg {
	static constexpr uint8_t fixed_exponent = 0; //set this to a value greater than 0 to have fixed size denominators
	
	//Maximum numerator size is 63 Bits resulting in 127 Bits for the sphere coordinate (which still fits a int128_t)
	struct Digits {
		static constexpr int NUMERATOR0=31;
		static constexpr int NUMERATOR1=NUMERATOR0;
		static constexpr int EXPONENT=7; //don't change this
// 		static constexpr int DENOMINATOR=static_cast<int>(1) << EXPONENT;
		static constexpr int DENOMINATOR=NUMERATOR0;
		static constexpr int POS=3; //don't change this
	};
	struct Signedness {
		static constexpr bool NUMERATOR0=true;
		static constexpr bool NUMERATOR1=true;
		static constexpr bool DENOMINATOR=false;
		static constexpr bool EXPONENT=DENOMINATOR;
		static constexpr bool POS=false;
	};
	//Bit sizes are not accurate for non integral types
	struct BitSizes {
		static constexpr int NUMERATOR0=Digits::NUMERATOR0+int(Signedness::NUMERATOR0);
		static constexpr int NUMERATOR1=Digits::NUMERATOR1+int(Signedness::NUMERATOR1);
		static constexpr int DENOMINATOR=Digits::DENOMINATOR+int(Signedness::DENOMINATOR);
		static constexpr int EXPONENT=Digits::EXPONENT+int(Signedness::EXPONENT);
		static constexpr int POS=Digits::POS+int(Signedness::POS);
	};
	
	struct Storage {
		static constexpr bool dense_coding = false;
		static constexpr int total_bits = (BitSizes::NUMERATOR0+BitSizes::NUMERATOR1+BitSizes::EXPONENT+BitSizes::POS);
		//if dense_coding==true then this tells the bit size, otherwise the byte size
		struct Size {
			static constexpr int NUMERATOR0 = dense_coding ? BitSizes::NUMERATOR0 : Kernel_sp::align(BitSizes::NUMERATOR0, 8);
			static constexpr int NUMERATOR1 = dense_coding ? BitSizes::NUMERATOR1 : Kernel_sp::align(BitSizes::NUMERATOR1, 8);
			static constexpr int POS = dense_coding ? BitSizes::POS : Kernel_sp::align(BitSizes::POS, 8);
			static constexpr int EXPONENT = fixed_exponent ? 0 : (dense_coding ? BitSizes::EXPONENT : Kernel_sp::align(BitSizes::EXPONENT, 8));
		};
		//if dense_coding==true then this tells the bit offset, otherwise the byte offset
		struct Begin {
			static constexpr int NUMERATOR0 = 0;
			static constexpr int NUMERATOR1 = NUMERATOR0+Size::NUMERATOR0;
			static constexpr int POS = NUMERATOR1+Size::NUMERATOR1;
			static constexpr int EXPONENT = POS+Size::POS;
		};
		//if dense_coding==true then this tells the bit offset, otherwise the byte offset
		struct End { //one passed the end
			static constexpr int NUMERATOR0 = Begin::NUMERATOR1;
			static constexpr int NUMERATOR1 = Begin::POS;
			static constexpr int POS = Begin::EXPONENT;
			static constexpr int EXPONENT = fixed_exponent > 0 ? POS+Size::EXPONENT : POS;
		};
		
		static constexpr int total_bytes = dense_coding ? Kernel_sp::align(total_bits, 8) : (Size::NUMERATOR0+Size::NUMERATOR1+Size::POS+Size::EXPONENT);
// 		static_assert(total_bytes == 10);
	};
	
	struct Types {
		using numerator = Kernel_sp::AINT<Digits::NUMERATOR0>;
		using denominator = Kernel_sp::AINT<Digits::DENOMINATOR>;
		
		using exponent = uint8_t;
		using pos = int8_t;
		
		using sphere_coord = Kernel_sp::AINT<2*Digits::NUMERATOR0+1>;
		using unsigned_sphere_coord = Kernel_sp::AUINT<2*Digits::NUMERATOR0+1>;
	};
	
	static constexpr Types::denominator fixed_denominator = (fixed_exponent > 0 ? (static_cast<Types::denominator>(1) << (fixed_exponent-1)) : Types::denominator(0));
	
	static_assert(std::is_integral<Types::numerator>::value, "Point_sp coordinates need to be integral values");
};

class Point_sp_base: public Point_sp_cfg {
public:
	using Config = Point_sp_cfg;
public:
	Point_sp_base();
	Point_sp_base(Point_sp_base const &) = default;
	Point_sp_base(Types::numerator _num0, Types::numerator _num1, Types::denominator _den, Types::pos _pos);
	virtual ~Point_sp_base();
	Point_sp_base & operator=(Point_sp_base const &) = default;
public:
	void set_numerator0(Types::numerator v);
	void set_numerator1(Types::numerator v);
	void set_denominator(Types::denominator v);
	void set_exponent(Types::exponent v);
	void set_pos(Types::pos v);
	void set_numerator0(mpz_class v);
	void set_numerator1(mpz_class v);
	void set_denominator(mpz_class v);
public:
	Types::numerator numerator0() const;
	Types::numerator numerator1() const;
	Types::denominator denominator() const;
	///@return position on Sphere
	LIB_RATSS_NAMESPACE::PositionOnSphere pos() const;
	Types::exponent exponent() const;
public:
	bool operator!=(Point_sp_base const & other) const;
	bool operator==(Point_sp_base const & other) const;
private:
	//Coordinates in the plane are of the form
	//x = num/2^s with s < 128
	//Thus we use 7 Bits to encode s
	//Additionally we need 3 Bits to encode the 6 position on the sphere
	//The remaining 84 Bits are used for the coordinates
	//Order is
	//numerator0, numerator1, pos, exponent
	std::array<char, Storage::total_bytes> m_d;
};

std::ostream & operator<<(std::ostream & out, Point_sp_base const & v);
	
} //end namespace detail

template<typename T_LINEAR_KERNEL>
class Point_sp: public detail::Point_sp_base {
public:
	using MyParent = detail::Point_sp_base;
	using Types = MyParent::Types;
	using MyBaseTrait = T_LINEAR_KERNEL;
	using FT = typename MyBaseTrait::FT;
	using Point_3 = typename MyBaseTrait::Point_3;
	using K = MyBaseTrait;
	using Numerators_3 = detail::Kernel_sp::Numerators_3<
		typename detail::Kernel_sp::AlignedIntegerTypeFromDigits<
			2*std::max<int>(
				fixed_exponent+1,
				MyParent::Config::Digits::NUMERATOR0
			)
		>::type
	>;
	using Numerator_Denominator_1_3 = detail::Kernel_sp::Numerator_Denominator_1_3<typename Numerators_3::RT>;
public:
	Point_sp();
	Point_sp(MyParent const & v);
	Point_sp(Types::numerator _num0, Types::numerator _num1, Types::denominator _den, Types::pos _pos);
	Point_sp(FT const & x, FT const & y, FT const & z);
	Point_sp(Point_3 const & v);
	Point_sp(Point_sp const & other);
	template<typename T_OTHER_TRAITS>
	Point_sp(Point_sp<T_OTHER_TRAITS> const & other);
public:
	Point_sp & operator=(Point_sp const & other) = default;
public:
	Point_3 point3() const;
	Numerators_3 numerators3() const;
	Numerator_Denominator_1_3 num_den_x_3() const;
	Numerator_Denominator_1_3 num_den_y_3() const;
	Numerator_Denominator_1_3 num_den_z_3() const;
public:
	template<typename T_FT = FT>
	T_FT x() const;
	
	template<typename T_FT = FT>
	T_FT y() const;
	
	template<typename T_FT = FT>
	T_FT z() const;
public:
	operator Point_3() const { return point3(); }
private:
	Point_3 point3_slow() const;
};

template<typename K1, typename K2>
bool operator==(Point_sp<K1> const & a, Point_sp<K2> const & b) {
	return static_cast<detail::Point_sp_base const &>(a) == static_cast<detail::Point_sp_base const &>(b);
}

template<typename K1, typename K2>
bool operator!=(Point_sp<K1> const & a, Point_sp<K2> const & b) {
	return static_cast<detail::Point_sp_base const &>(a) != static_cast<detail::Point_sp_base const &>(b);
}

template<typename T_BASE_TRAITS>
std::ostream & operator<<(std::ostream & out, Point_sp<T_BASE_TRAITS> const & v) {
	out << static_cast<detail::Point_sp_base const &>(v);
	out << "=(" << v.x() << ", " << v.y() << ", " << v.z() << ")";
	return out;
}

}//end namespace

//the definitions

namespace LIB_DTS2_NAMESPACE {

#define PTSP_TMP_PRMS template<typename T_BASE_TRAITS>
#define PTSP_CLS_NAME Point_sp<T_BASE_TRAITS>

PTSP_TMP_PRMS
PTSP_CLS_NAME::Point_sp() {}

PTSP_TMP_PRMS
PTSP_CLS_NAME::Point_sp(Point_sp const & other) :
MyParent(other)
{}

PTSP_TMP_PRMS
template<typename T_OTHER_TRAITS>
PTSP_CLS_NAME::Point_sp(Point_sp<T_OTHER_TRAITS> const & other) :
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
	
	if (MyParent::fixed_exponent) {
		if (mden > MyParent::fixed_denominator) {
			throw std::runtime_error("Point_sp: Input coordinate precision too high");
		}
		else {
			mden = MyParent::fixed_denominator;
		}
	}
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
	assert(this->x() == x);
	assert(this->y() == y);
	assert(this->z() == z);
	assert(Point_3(x,y,z) == point3());
}

PTSP_TMP_PRMS
PTSP_CLS_NAME::Point_sp(Types::numerator _num0, Types::numerator _num1, Types::denominator _den, Types::pos _pos) :
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
	if (exponent() > 31) {
		return point3_slow();
	}
	using LIB_RATSS_NAMESPACE::convert;
// 	using PQType = mpq_class;
	auto PQType = [](Types::sphere_coord const & x) -> FT { return LIB_RATSS_NAMESPACE::convert<FT>(x); };
	FT x, y, z;
	auto sqr = [](Types::sphere_coord v) { return v*v; };
	Types::sphere_coord den2 = sqr(denominator());
	Types::sphere_coord sum_p_i2 = sqr(numerator0()) + sqr(numerator1());
	switch (abs(pos())) {
	case 1: //x is the missing coordinate, thus num0 holds y, num1 holds z
		x = (std::signbit<int>(pos()) ? 1 : -1)*PQType(sum_p_i2 - den2)/PQType(den2 + sum_p_i2);
		y = PQType(2*numerator0()*denominator())/PQType(den2 + sum_p_i2);
		z = PQType(2*numerator1()*denominator())/PQType(den2 + sum_p_i2);
		break;
	case 2: //y is the missing coordinate, thus num0 holds x, num1 holds z
		x = PQType(2*numerator0()*denominator())/PQType(den2 + sum_p_i2);
		y = (std::signbit<int>(pos()) ? 1 : -1)*PQType(sum_p_i2 - den2)/PQType(den2 + sum_p_i2);
		z = PQType(2*numerator1()*denominator())/PQType(den2 + sum_p_i2);
		break;
	case 3: //z is the missing coordinate, thus num0 holds x, num1 holds y
		x = PQType(2*numerator0()*denominator())/PQType(den2 + sum_p_i2);
		y = PQType(2*numerator1()*denominator())/PQType(den2 + sum_p_i2);
		z = (std::signbit<int>(pos()) ? 1 : -1)*PQType(sum_p_i2 - den2)/PQType(den2 + sum_p_i2);
		break;
	};
	Point_3 result(convert<FT>(x), convert<FT>(y), convert<FT>(z));
	assert(result == point3_slow());
	return result;
}

PTSP_TMP_PRMS
typename PTSP_CLS_NAME::Numerators_3
PTSP_CLS_NAME::numerators3() const {
	if (pos() == 0) {  //Hack for LIB_DTS2_ORIGIN 
		return Numerators_3(LIB_DTS2_ORIGIN_X, LIB_DTS2_ORIGIN_Y, LIB_DTS2_ORIGIN_Z);
	}
	using RT = typename Numerators_3::RT;
	if (2*exponent()+1 > std::numeric_limits<RT>::digits) {
		throw std::runtime_error("Point_sp::numerators3: exponent is too large");
	}
	using LIB_RATSS_NAMESPACE::convert;
	Numerators_3 r;
	auto sqr = [](Types::sphere_coord v) { return v*v; };
	Types::sphere_coord den2 = sqr(denominator());
	Types::sphere_coord sum_p_i2 = sqr(numerator0()) + sqr(numerator1());
	switch (abs(pos())) {
	case 1: //x is the missing coordinate, thus num0 holds y, num1 holds z
		r.x = (std::signbit<int>(pos()) ? 1 : -1)*RT(sum_p_i2 - den2);
		r.y = RT(2*numerator0()*denominator());
		r.z = RT(2*numerator1()*denominator());
		break;
	case 2: //y is the missing coordinate, thus num0 holds x, num1 holds z
		r.x = RT(2*numerator0()*denominator());
		r.y = (std::signbit<int>(pos()) ? 1 : -1)*RT(sum_p_i2 - den2);
		r.z = RT(2*numerator1()*denominator());
		break;
	case 3: //z is the missing coordinate, thus num0 holds x, num1 holds y
		r.x = RT(2*numerator0()*denominator());
		r.y = RT(2*numerator1()*denominator());
		r.z = (std::signbit<int>(pos()) ? 1 : -1)*RT(sum_p_i2 - den2);
		break;
	};
	return r;
}

PTSP_TMP_PRMS
typename PTSP_CLS_NAME::Numerator_Denominator_1_3
PTSP_CLS_NAME::num_den_x_3() const {
	if (pos() == 0) {  //Hack for LIB_DTS2_ORIGIN 
		return Numerator_Denominator_1_3(LIB_DTS2_ORIGIN_X, 1);
	}
	using RT = typename Numerator_Denominator_1_3::RT;
	if (2*exponent()+1 > std::numeric_limits<RT>::digits) {
		throw std::runtime_error("Point_sp::num_den_x_3: exponent is too large");
	}
	auto sqr = [](Types::sphere_coord v) { return v*v; };
	Types::sphere_coord den2 = sqr(denominator());
	Types::sphere_coord sum_p_i2 = sqr(numerator0()) + sqr(numerator1());
	Numerator_Denominator_1_3 result;
	switch (abs(pos())) {
	case 1: //x is the missing coordinate
		result = Numerator_Denominator_1_3((std::signbit<int>(pos()) ? 1 : -1)*RT(sum_p_i2 - den2), RT(den2 + sum_p_i2));
		break;
	default: //in any case num0 holds x
		result = Numerator_Denominator_1_3(RT(2*numerator0()*denominator()), RT(den2 + sum_p_i2));
		break;
	};
	{
		using LIB_RATSS_NAMESPACE::convert;
		assert(convert<mpq_class>(result.num)/convert<mpq_class>(result.den) == x<mpq_class>());
	}
	return result;
}

PTSP_TMP_PRMS
typename PTSP_CLS_NAME::Numerator_Denominator_1_3
PTSP_CLS_NAME::num_den_y_3() const {
	if (pos() == 0) {  //Hack for LIB_DTS2_ORIGIN 
		return Numerator_Denominator_1_3(LIB_DTS2_ORIGIN_Y, 1);
	}
	using RT = typename Numerator_Denominator_1_3::RT;
	if (2*exponent()+1 > std::numeric_limits<RT>::digits) {
		throw std::runtime_error("Point_sp::num_den_y_3: exponent is too large");
	}
	auto sqr = [](Types::sphere_coord v) { return v*v; };
	Types::sphere_coord den2 = sqr(denominator());
	Types::sphere_coord sum_p_i2 = sqr(numerator0()) + sqr(numerator1());
	Numerator_Denominator_1_3 result;
	switch (abs(pos())) {
	case 1: //x is the missing coordinate, thus num0 holds y
		result = Numerator_Denominator_1_3(RT(2*numerator0()*denominator()), RT(den2 + sum_p_i2));
		break;
	case 2: //y is the missing coordinate
		result = Numerator_Denominator_1_3((std::signbit<int>(pos()) ? 1 : -1)*RT(sum_p_i2 - den2), RT(den2 + sum_p_i2));
		break;
	case 3: //z is the missing coordinate, thus num1 holds y
		result = Numerator_Denominator_1_3(RT(2*numerator1()*denominator()), RT(den2 + sum_p_i2));
		break;
	};
	{
		using LIB_RATSS_NAMESPACE::convert;
		assert(convert<mpq_class>(result.num)/convert<mpq_class>(result.den) == y<mpq_class>());
	}
	return result;
}

PTSP_TMP_PRMS
typename PTSP_CLS_NAME::Numerator_Denominator_1_3
PTSP_CLS_NAME::num_den_z_3() const {
	if (pos() == 0) {  //Hack for LIB_DTS2_ORIGIN 
		return Numerator_Denominator_1_3(LIB_DTS2_ORIGIN_Z, 1);
	}
	using RT = typename Numerator_Denominator_1_3::RT;
	if (2*exponent()+1 > std::numeric_limits<RT>::digits) {
		throw std::runtime_error("Point_sp::num_den_z_3: exponent is too large");
	}
	auto sqr = [](Types::sphere_coord v) { return v*v; };
	Types::sphere_coord den2 = sqr(denominator());
	Types::sphere_coord sum_p_i2 = sqr(numerator0()) + sqr(numerator1());
	Numerator_Denominator_1_3 result;
	switch (abs(pos())) {
	case 1: //x,y are the missing coordinate, thus num1 holds z
	case 2:
		result = Numerator_Denominator_1_3(RT(2*numerator1()*denominator()), RT(den2 + sum_p_i2));
		break;
	case 3: //z is the missing coordinate, thus num1 holds y
		result = Numerator_Denominator_1_3((std::signbit<int>(pos()) ? 1 : -1)*RT(sum_p_i2 - den2), RT(den2 + sum_p_i2));
		break;
	};
	{
		using LIB_RATSS_NAMESPACE::convert;
		assert(convert<mpq_class>(result.num)/convert<mpq_class>(result.den) == z<mpq_class>());
	}
	return result;
}

PTSP_TMP_PRMS
typename PTSP_CLS_NAME::Point_3
PTSP_CLS_NAME::point3_slow() const {
	
	LIB_RATSS_NAMESPACE::ProjectS2 proj;
	mpq_class xp, yp, zp;
	mpq_class xs, ys, zs;
	mpq_class num0 = LIB_RATSS_NAMESPACE::convert<mpq_class>(numerator0());
	mpq_class num1 = LIB_RATSS_NAMESPACE::convert<mpq_class>(numerator1());
	mpq_class den = LIB_RATSS_NAMESPACE::convert<mpq_class>(denominator());
	
	switch (std::abs(this->pos())) {
	case 1:
		xp = 0;
		yp = num0/den;
		zp = num1/den;
		break;
	case 2:
		xp = num0/den;
		yp = 0;
		zp = num1/den;
		break;
	case 3:
		xp = num0/den;
		yp = num1/den;
		zp = 0;
		break;
	default:
		throw std::runtime_error("Point_sp: Position is out of bounds");
	};
	proj.plane2Sphere(xp, yp, zp, pos(), xs, ys, zs);
	
	using LIB_RATSS_NAMESPACE::convert;
	return Point_3( convert<FT>(xs), convert<FT>(ys), convert<FT>(zs) );
}

//plane2sphere is as follows:
// 2x_i/(1 + sum x_i^2)
//=2(p_i/den)/(1+1/den^2 * sum p_i^2)
//=2(p_i/den)/(1/den^2*(den^2 + sum p_i^2)
//=2*p_i*den/(den^2 + sum p_i^2)
//
//Or for the removed dimension
// (sum x_i^2 - 1)/(1 + sum x_i^2)
//=(1/den^2*(sum p_i^2 - den^2))/(1/den^2*(den^2 + sum p_i^2))
//=(sum p_i^2 - den^2)/(den^2 + sum p_i^2)
//We have that {den, p_i} <= 2^31
//Thus den^2 <= 2^62 and den^2 + den^2 <= 2^63
//Hence int64_t should be enough


PTSP_TMP_PRMS
template<typename T_FT>
T_FT
PTSP_CLS_NAME::x() const {
	if (pos() == 0) {  //Hack for LIB_DTS2_ORIGIN 
		return T_FT(LIB_DTS2_ORIGIN_X);
	}
	if (exponent() > 31) {
		return ratss::convert<T_FT>( point3().x() );
	}
	using LIB_RATSS_NAMESPACE::convert;
	auto sqr = [](Types::sphere_coord v) { return v*v; };
	Types::sphere_coord den2 = sqr(denominator());
	Types::sphere_coord sum_p_i2 = sqr(numerator0()) + sqr(numerator1());
	mpq_class result;
	switch (abs(pos())) {
	case 1: //x is the missing coordinate
		result = (std::signbit<int>(pos()) ? 1 : -1)*convert<mpq_class>(sum_p_i2 - den2)/convert<mpq_class>(den2 + sum_p_i2);
		break;
	default: //in any case num0 holds x
		result = convert<mpq_class>(2*numerator0()*denominator())/convert<mpq_class>(den2 + sum_p_i2);
		break;
	};
	assert(result == convert<mpq_class>(point3_slow().x()));
	return convert<T_FT>(result);
}

PTSP_TMP_PRMS
template<typename T_FT>
T_FT
PTSP_CLS_NAME::y() const {
	if (pos() == 0) {  //Hack for LIB_DTS2_ORIGIN 
		return T_FT(LIB_DTS2_ORIGIN_Y);
	}
	if (exponent() > 31) {
		return ratss::convert<T_FT>( point3().y() );
	}
	using LIB_RATSS_NAMESPACE::convert;
	auto sqr = [](Types::sphere_coord v) { return v*v; };
	Types::sphere_coord den2 = sqr(denominator());
	Types::sphere_coord sum_p_i2 = sqr(numerator0()) + sqr(numerator1());
	mpq_class result;
	switch (abs(pos())) {
	case 1: //x is the missing coordinate, thus num0 holds y
		result = convert<mpq_class>(2*numerator0()*denominator())/convert<mpq_class>(den2 + sum_p_i2);
		break;
	case 2: //y is the missing coordinate
		result = (std::signbit<int>(pos()) ? 1 : -1)*convert<mpq_class>(sum_p_i2 - den2)/convert<mpq_class>(den2 + sum_p_i2);
		break;
	case 3: //z is the missing coordinate, thus num1 holds y
		result = convert<mpq_class>(2*numerator1()*denominator())/convert<mpq_class>(den2 + sum_p_i2);
		break;
	};
	assert(result == convert<mpq_class>(point3_slow().y()));
	return convert<T_FT>(result);
}

PTSP_TMP_PRMS
template<typename T_FT>
T_FT
PTSP_CLS_NAME::z() const {
	if (pos() == 0) {  //Hack for LIB_DTS2_ORIGIN 
		return T_FT(LIB_DTS2_ORIGIN_Z);
	}
	if (exponent() > 31) {
		return ratss::convert<T_FT>( point3().z() );
	}
	using LIB_RATSS_NAMESPACE::convert;
	auto sqr = [](Types::sphere_coord v) { return v*v; };
	Types::sphere_coord den2 = sqr(denominator());
	Types::sphere_coord sum_p_i2 = sqr(numerator0()) + sqr(numerator1());
	mpq_class result;
	switch (abs(pos())) {
	case 1: //x,y are the missing coordinate, thus num1 holds z
	case 2:
		result = convert<mpq_class>(2*numerator1()*denominator())/convert<mpq_class>(den2 + sum_p_i2);
		break;
	case 3: //z is the missing coordinate, thus num1 holds y
		result = (std::signbit<int>(pos()) ? 1 : -1)*convert<mpq_class>(sum_p_i2 - den2)/convert<mpq_class>(den2 + sum_p_i2);
		break;
	};
	assert(result == convert<mpq_class>(point3_slow().z()));
	return convert<T_FT>(result);
}

#undef PTSP_TMP_PRMS
#undef PTSP_CLS_NAME

};

#endif
