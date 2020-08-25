#pragma once
#ifndef LIB_DTS2_CONSTRAINED_DELAUNAY_TRIANGULATION_BASE_TRAITS_S2_H
#define LIB_DTS2_CONSTRAINED_DELAUNAY_TRIANGULATION_BASE_TRAITS_S2_H

#include <libdts2/constants.h>
#include <libdts2/util.h>
#include <libratss/ProjectS2.h>
#include <libratss/SphericalCoord.h>
#include <libratss/GeoCoord.h>

#include <CGAL/number_utils.h>

#include <CGAL/enum.h>


#include <assert.h>

namespace LIB_DTS2_NAMESPACE {
namespace detail {

template<typename T_LINEAR_KERNEL>
class EpsBasedAuxPoints {
public:
	using BaseTraits = T_LINEAR_KERNEL;
	using Point_3 = typename BaseTraits::Point_3;
	using FT = typename BaseTraits::FT;
public:
	class Generate_auxiliary_point final {
	public:
		Generate_auxiliary_point(FT const & _epsZ) : m_epsZ(_epsZ) {}
	public:
		inline Point_3 operator()(AuxPointSelector s) const {
			if (s == AuxPointSelector::LOWER) {
				return Point_3(0, 0, -1);
			}
			FT epsilon(LIB_RATSS_NAMESPACE::Conversion<FT>::moveFrom(mpq_class(1, 0x8000)));
			switch(s) {
			case AuxPointSelector::UPPER_0:
				return generateAuxPoint(0, epsilon);
			case AuxPointSelector::UPPER_1:
				return generateAuxPoint(-epsilon, -epsilon/2);
			case AuxPointSelector::UPPER_2:
				return generateAuxPoint(epsilon, -epsilon/2);
			default:
				throw std::runtime_error("Invalid AuxPointSelector value");
			};
		};
	protected:
		Point_3 generateAuxPoint(FT xp, FT yp) const {
			FT dummy(0);
			FT xs, ys, zs, zs_prev;
			m_p.plane2Sphere(xp, yp, dummy, LIB_RATSS_NAMESPACE::SP_UPPER, xs, ys, zs);
			zs_prev = zs;
			while(zs < m_epsZ) {
				xp /= 2;
				yp /= 2;
				zs_prev = zs;
				m_p.plane2Sphere(xp, yp, dummy, LIB_RATSS_NAMESPACE::SP_UPPER, xs, ys, zs);
				assert(zs_prev < zs);
				assert(zs < 1);
			}
			return Point_3(xs, ys, zs);
		}
	private:
		LIB_RATSS_NAMESPACE::ProjectS2 m_p;
		FT m_epsZ;
	};
	
	class Is_auxiliary_point final {
	public:
		Is_auxiliary_point(FT const & _epsZ) : m_epsZ(_epsZ) {}
	public:
		template<typename T_POINT>
		inline bool operator()(T_POINT const & p) const {
			return p.z() == FT(-1) || p.z() >= m_epsZ;
		};
	private:
		FT m_epsZ;
	};
	
	class Is_in_auxiliary_triangle final {
	public:
		Is_in_auxiliary_triangle(FT const & _epsZ) : m_epsZ(_epsZ) {}
	public:
		template<typename T_POINT>
		inline bool operator()(T_POINT const & p) const {
			return p.z() >= m_epsZ;
		};
	private:
		FT m_epsZ;
	};
	
	//Is true if the given point is on the sphere and outside of the auxiliary triangle
	class Is_valid_point_on_sphere final {
	public:
		Is_valid_point_on_sphere(Is_auxiliary_point const & iap) : m_iap(iap) {}
	public:
		bool operator()(Point_3 const & p) const {
			return (p.x()*p.x()+p.y()*p.y()+p.z()*p.z()) == 1 && !m_iap(p);
		}
	private:
		Is_auxiliary_point m_iap;
	};
public:
	EpsBasedAuxPoints() :
		m_epsZ(FT(1)-FT(std::numeric_limits<double>::epsilon()))
	{}
	EpsBasedAuxPoints(FT const & _epsZ) :
		m_epsZ(_epsZ)
	{}
	virtual ~EpsBasedAuxPoints() {}
public:
	Is_auxiliary_point is_auxiliary_point_object() const {
		return Is_auxiliary_point(epsZ());
	}
	
	Is_in_auxiliary_triangle is_in_auxiliary_triangle_object() const {
		return Is_in_auxiliary_triangle(epsZ());
	}
	
	Is_valid_point_on_sphere is_valid_point_on_sphere_object() const {
		return Is_valid_point_on_sphere(is_auxiliary_point_object());
	}
	
	Generate_auxiliary_point generate_auxiliary_point_object() const {
		return Generate_auxiliary_point(epsZ());
	}
public:
	FT const & epsZ() const { return m_epsZ; }
private:
	FT m_epsZ;
};

//User defines 3 Points as the upper aux points and one point as the lower aux point
//Let ot=Orientation_3(u0, u1, u2, l) be the orientation of the lower point with respect to the upper points
//Then all points p with with Orientation_3(u0, u1, u2, p) != ot are within in the auxiliary triangle and
//thus outside of the resulting triangulation
template<typename T_LINEAR_KERNEL>
class UserDefinedAuxPoints {
public:
	using BaseTraits = T_LINEAR_KERNEL;
	using Point_3 = typename BaseTraits::Point_3;
	using FT = typename BaseTraits::FT;
private:
	struct Data {
		std::array<Point_3, 4> ap;
	};
public:
	class Generate_auxiliary_point final {
	public:
		Generate_auxiliary_point(std::shared_ptr<Data> _d) : m_d(_d) {}
	public:
		inline Point_3 operator()(AuxPointSelector s) const {
			return m_d->ap.at(static_cast<int>(s));
		};
	private:
		std::shared_ptr<Data> m_d;
	};
	
	class Is_auxiliary_point final {
	public:
		Is_auxiliary_point(std::shared_ptr<Data> _d) : m_d(_d) {}
	public:
		inline bool operator()(Point_3 const & p) const {
			for(auto const & x : m_d->ap) {
				if (x == p) {
					return true;
				}
			}
			return false;
		};
	private:
		typename BaseTraits::Orientation_3 m_ot3;
		std::shared_ptr<Data> m_d;
	};
	
	class Is_in_auxiliary_triangle final {
	public:
		Is_in_auxiliary_triangle(std::shared_ptr<Data> _d) : m_d(_d) {}
	public:
		inline bool operator()(Point_3 const & p) const {
			auto otp = m_ot3(m_d->ap[0], m_d->ap[1], m_d->ap[2], p);
			if (otp == CGAL::COLLINEAR) {
				return true;
			}
			return otp != m_ot3(m_d->ap[0], m_d->ap[1], m_d->ap[2], m_d->ap[3]);
		};
	private:
		typename BaseTraits::Orientation_3 m_ot3;
		std::shared_ptr<Data> m_d;
	};
	
	//Is true if the given point is on the sphere and outside of the auxiliary triangle
	class Is_valid_point_on_sphere final {
	public:
		Is_valid_point_on_sphere(Is_auxiliary_point const & iap) : m_iap(iap) {}
	public:
		bool operator()(Point_3 const & p) const {
			return (p.x()*p.x()+p.y()*p.y()+p.z()*p.z()) == 1 && !m_iap(p);
		}
	private:
		Is_auxiliary_point m_iap;
	};
public:
	UserDefinedAuxPoints() {}
	UserDefinedAuxPoints(Point_3 const & u0, Point_3 const & u1, Point_3 const & u2, Point_3 const & l) :
	m_d(std::make_shared<Data>())
	{
		m_d->ap[0] = u0;
		m_d->ap[1] = u1;
		m_d->ap[2] = u2;
		m_d->ap[3] = l;
	}
	virtual ~UserDefinedAuxPoints() {}
public:
	Is_auxiliary_point is_auxiliary_point_object() const {
		return Is_auxiliary_point(m_d);
	}
	
	Is_in_auxiliary_triangle is_in_auxiliary_triangle_object() const {
		return Is_in_auxiliary_triangle(m_d);
	}
	
	Is_valid_point_on_sphere is_valid_point_on_sphere_object() const {
		return Is_valid_point_on_sphere(is_auxiliary_point_object());
	}
	
	Generate_auxiliary_point generate_auxiliary_point_object() const {
		return Generate_auxiliary_point(m_d);
	}
private:
	std::shared_ptr<Data> m_d;
};

} //end namespace detail

#ifndef NDEBUG
	#define DEBUG_OUT(__STR) {std::cout <<  __STR << std::endl;}
#else
	#define DEBUG_OUT(__) {}
#endif

template<
	typename T_LINEAR_KERNEL,
	typename T_AUX_POINT_GENERATOR=detail::EpsBasedAuxPoints<T_LINEAR_KERNEL>,
	typename T_POINT = typename T_LINEAR_KERNEL::Point_3
>
class Constrained_delaunay_triangulation_base_traits_s2 {
public:
	using LinearKernel = T_LINEAR_KERNEL;
	using MyBaseTrait = T_LINEAR_KERNEL;
	using MyKernel = T_LINEAR_KERNEL;
	using MySelf = Constrained_delaunay_triangulation_base_traits_s2<MyBaseTrait>;
	using AuxiliaryPointsGenerator = T_AUX_POINT_GENERATOR;
protected: //we only want to expose typedefs that are really needed by the triangulation
	using Segment_3 = typename MyBaseTrait::Segment_3;
	using Orientation_3 = typename MyBaseTrait::Orientation_3;

	using Compare_distance_3 = typename MyBaseTrait::Compare_distance_3;
	
	//needed by Compare_x_2 and Compare_y_2
	using Compare_x_3 = typename MyBaseTrait::Compare_x_3;
	using Compare_y_3 = typename MyBaseTrait::Compare_y_3;
	using Compare_z_3 = typename MyBaseTrait::Compare_z_3;
	
	
	using Triangle_3 = typename MyBaseTrait::Triangle_3;
	
	//stuff for intersection support
	using Construct_segment_3 = typename MyBaseTrait::Construct_segment_3;
	using Construct_point_3 = typename MyBaseTrait::Construct_point_3;
	using Intersect_3 = typename MyBaseTrait::Intersect_3;
	using Do_intersect_3 = typename MyBaseTrait::Do_intersect_3;
	using Plane_3 = typename MyBaseTrait::Plane_3;
	using Line_3 = typename MyBaseTrait::Line_3;
	using Ray_3 = typename MyBaseTrait::Ray_3;
	using Vector_3 = typename MyBaseTrait::Vector_3;
	
	//stuff used for Collinear_are_ordered_along_line_2
	using Coplanar_3 = typename MyBaseTrait::Coplanar_3;
	using Compute_scalar_product_3 = typename MyBaseTrait::Compute_scalar_product_3;
	using Construct_cross_product_vector_3 = typename MyBaseTrait::Construct_cross_product_vector_3;
	
public:
	using FT = typename MyBaseTrait::FT;
	using Oriented_side = typename MyBaseTrait::Oriented_side;
	using Comparison_result = typename MyBaseTrait::Comparison_result;
	using Orientation = typename MyBaseTrait::Orientation;
	
	using Segment_2 = Segment_3;
	using Segment = Segment_2;
	
	//we want our triangulation to live in cartesian 3D-space
	//the CDT lives in 2D euclidean space so we have to redefine Point_2;
	using Point_3 = T_POINT; //typename MyBaseTrait::Point_3;
	using Point_2 = T_POINT;
	using Point = Point_2;
	
	using Compare_distance_2 = Compare_distance_3;
	
	//the following is needed for spatial sorting
	using Less_x_3 = typename MyBaseTrait::Less_x_3;
	using Less_y_3 = typename MyBaseTrait::Less_y_3;
	using Less_z_3 = typename MyBaseTrait::Less_z_3;
	
	//this is needed for spatial sorting during insertion of constraints
	//cgal uses Spatial_sort_traits_adapter_2 to sort inserted vertices
	//this in turn needs the following predicates
	//maybe we should just specialize Spatial_sort_traits_adapter_2 for this Trait
	
	using Less_x_2 = Less_x_3;
	using Less_y_2 = Less_y_3;
	
	using Construct_segment_2 = Construct_segment_3;
	using Construct_point_2 = Construct_point_3;
public:
	using Projector = LIB_RATSS_NAMESPACE::ProjectS2;
public:
	using Is_auxiliary_point = typename AuxiliaryPointsGenerator::Is_auxiliary_point;
	using Is_in_auxiliary_triangle = typename AuxiliaryPointsGenerator::Is_in_auxiliary_triangle;
	using Is_valid_point_on_sphere = typename AuxiliaryPointsGenerator::Is_valid_point_on_sphere;
	using Generate_auxiliary_point = typename AuxiliaryPointsGenerator::Generate_auxiliary_point;
protected:
	template<typename T>
	static auto bpc(T && x) {
		using BasePoint = typename LinearKernel::Point_3;
		return detail::ipt_static_cast<BasePoint>(std::forward<T>(x));
	}
public: //own implementations

	class Triangle_2: private Triangle_3 {
	public:
		Triangle_2() {}
		Triangle_2(const Point_3 &p, const Point_3 &q, const Point_3 &r) : Triangle_3(p, q, r) {
			assert(false);
		}
		bool has_on(const Point_3 &p) const {
			assert(false);
			return Triangle_3::has_on(p);
		}
	};

	class Orientation_2
	{
	public:
		Orientation_2(Orientation_3 const & _ot3, Is_in_auxiliary_triangle const & _iat) : ot3(_ot3), m_iat(_iat) {}
		Orientation operator()(const Point_3 & p, const Point_3 & q, const Point_3 & r) const
		{
			assert(p.z() <= 1 && q.z() <= 1 && r.z() <= 1);
			assert(p != q);
			
			Orientation oriented_side = ot3(bpc(p), bpc(q), bpc(LIB_DTS2_ORIGIN), bpc(r));
			
			if (oriented_side != Orientation::COLLINEAR && m_iat(p) && m_iat(q)) {
				//the point always has to be on the opposite side of the infinite vertex
				oriented_side = - ot3(bpc(p), bpc(q), bpc(LIB_DTS2_ORIGIN), bpc(Point(0, 0, 1)));
			}
			return oriented_side;
		}
	private:
		Orientation_3 ot3;
		Is_in_auxiliary_triangle m_iat;
	};
	
	class Side_of_oriented_circle_2 {
	public:
		Side_of_oriented_circle_2(const Orientation_3 & _ot3) : ot3(_ot3) {}
		Orientation operator()(const Point_3 & p, const Point_3 & q, const Point_3 & r, const Point_3 & s) const
		{
			assert(p != q && p != r && q != r && r != s);
			return ot3(bpc(p), bpc(q), bpc(r), bpc(s));
		}
	private:
		Orientation_3 ot3;
	};
	
	class Compare_x_2 {
	public:
		Compare_x_2() {}
	public:
		Comparison_result  operator()(const Point_3 & p, const Point_3 & q) const {
			return m_cx3(bpc(p), bpc(q));
		}
	private:
		Compare_x_3 m_cx3;
		Compare_z_3 m_cz3;
	};
	
	class Compare_y_2 {
	public:
		Compare_y_2() {}
	public:
		Comparison_result  operator()(const Point_3 & p, const Point_3 & q) const {
			return m_cy3(bpc(p), bpc(q));
		}
	private:
		Compare_y_3 m_cy3;
	};
	
	///Checks if 3 points p, q,r are ordered along a great circle
	///They are ordered along the great circle if the angle >pq is smaller than >pr
	///AND if r,q are on the same side of the line from p to -p
	class Collinear_are_ordered_along_line_2 {
	public:
		Collinear_are_ordered_along_line_2() {}
	public:
		bool operator()(const Point_3 & p, const Point_3 & q, const Point_3 & r) const {
			if (!m_cop3(bpc(LIB_DTS2_ORIGIN), bpc(p), bpc(q), bpc(r))) {
				return false;
			}
			if (p == q || q == r) {
				return true;
			}
			//are coplanar, the direction of one point to another is defined by the vector with shorter euclidean distance
			//the following works since p, q, r have all the same length
			//here we check the angle >pq and >pr. >pq has to be smaller than >pr
			auto pq_sp = m_csp3(bpc(p), bpc(q));
			auto pr_sp = m_csp3(bpc(p), bpc(r));
			if (pq_sp > pr_sp) {
				return false;
			}
			//>pr has a larger angle than >pq BUT it can still be on the other other side of the circle
			//we have to make sure that r, q are on the same side of the line from p to -p
			
			//if the cross-product of pq and pr point the same direction,
			//then q,r are on the same side
			//They point in the same direction if pq_cpv = l*pr_cpv with l positive
			//They point in different direction if l negative
			auto pq_cpv = m_ccpv3(bpc(p), bpc(q));
			auto pr_cpv = m_ccpv3(bpc(p), bpc(r));
			return (CGAL::sign(pq_cpv.x()) * CGAL::sign(pr_cpv.x()) == CGAL::POSITIVE ||
					CGAL::sign(pq_cpv.y()) * CGAL::sign(pr_cpv.y()) == CGAL::POSITIVE ||
					CGAL::sign(pq_cpv.z()) * CGAL::sign(pr_cpv.z()) == CGAL::POSITIVE);
		}
	private:
		Coplanar_3 m_cop3;
		Compute_scalar_product_3 m_csp3;
		Construct_cross_product_vector_3 m_ccpv3;
	};
	
public:
	
	///sub classes of this trait should sub-class this class and add an appropriate
	///Point_3 operator()(const Point_3 & p) const; function to map points to the sphere 
	class Project_on_sphere {
	private:
		template<typename T>
		using Conversion = LIB_RATSS_NAMESPACE::Conversion<T>;
	public:
		Project_on_sphere(int _significands, const Projector & p) : m_significands(_significands), m_proj(p) {}
		Project_on_sphere(int _significands, int _snapType, const Projector & p) : m_significands(_significands), m_snapType(_snapType), m_proj(p) {}
	public:
		Point_3 operator()(const std::pair<double, double> & geoCoord) const {
			return this->operator()( LIB_RATSS_NAMESPACE::GeoCoord(geoCoord) );
		}
		Point_3 operator()(const LIB_RATSS_NAMESPACE::GeoCoord & geoCoord) const {
			mpq_class x, y, z;
			m_proj.projectFromGeo(mpfr::mpreal(geoCoord.lat), mpfr::mpreal(geoCoord.lon), x, y, z, significands(), snapType());
			return Point_3(
				Conversion<FT>::moveFrom(x),
				Conversion<FT>::moveFrom(y),
				Conversion<FT>::moveFrom(z)
			);
		}
		Point_3 operator()(const LIB_RATSS_NAMESPACE::SphericalCoord & sphericalCoord) const {
			mpq_class x, y, z;
			m_proj.projectFromSpherical(mpfr::mpreal(sphericalCoord.theta), mpfr::mpreal(sphericalCoord.phi), x, y, z, significands(), snapType());
			return Point_3(
				Conversion<FT>::moveFrom(x),
				Conversion<FT>::moveFrom(y),
				Conversion<FT>::moveFrom(z)
			);
		}
	public:
		int significands() const {
			return m_significands;
		}
		int snapType() const {
			return m_snapType;
		}
		const Projector & projector() const {
			return m_proj;
		}
	private:
		int m_significands;
		int m_snapType{LIB_RATSS_NAMESPACE::ST_FX|LIB_RATSS_NAMESPACE::ST_NORMALIZE|LIB_RATSS_NAMESPACE::ST_PLANE};
		Projector m_proj;
	};
protected:
	
	Orientation_3 orientation_3_object() const {
		return baseTraits().orientation_3_object();
	}
	
	//needed for intersections
	
	Construct_segment_3 construct_segment_3_object() const {
		return baseTraits().construct_segment_3_object();
	}

	Construct_point_3 construct_point_3_object() const {
		return baseTraits().construct_point_3_object();
	}
	
	Intersect_3 intersect_3_object() const {
		return baseTraits().intersect_3_object();
	}
	
	Do_intersect_3 do_intersect_3_object() const {
		return baseTraits().do_intersect_3_object();
	}
protected:
	static void print(const Point_3 & p) {
		using Conversion = LIB_RATSS_NAMESPACE::Conversion<FT>;
		std::cout << '(' << doubleValue(p.x()) << ','
		<< doubleValue(p.y()) << ','
		<< doubleValue(p.z()) << ")=("
		<< Conversion::toMpq(p.x()) << ", " << Conversion::toMpq(p.y()) << ", " << Conversion::toMpq(p.z()) << ")\n";
	}
	static void print(std::ostream & out, const Point_3 & p) {
		using Conversion = LIB_RATSS_NAMESPACE::Conversion<FT>;
		out << '(' << Conversion::toMpq(p.x()) << ", " << Conversion::toMpq(p.y()) << ", " << Conversion::toMpq(p.z()) << ')';
	}
public:
	const FT & epsilon() const;
public:
	static double doubleValue(const FT & v) {
		return LIB_RATSS_NAMESPACE::Conversion<FT>::toMpreal(v, 53).toDouble();
	}
	static mpq_class mpqValue(const FT & v) {
		return LIB_RATSS_NAMESPACE::Conversion<FT>::toMpq(v);
	}
	static std::string ratString(const FT & v) {
		std::stringstream ss;
		ss << mpqValue(v);
		return ss.str();
	}
public: //object functions
	///This does not correctly initialize this trait!
	Constrained_delaunay_triangulation_base_traits_s2() {}
	Constrained_delaunay_triangulation_base_traits_s2(int _significands) :
	m_significands(_significands)
	{
		assert(m_significands > 2);
	}
	///@param epsilon set the value of the z-coordinate above which no points should exist
	Constrained_delaunay_triangulation_base_traits_s2(AuxiliaryPointsGenerator const & _apg, int _significands) :
	m_apg(_apg),
	m_significands(_significands)
	{
		assert(m_significands > 2);
	}
	
	Constrained_delaunay_triangulation_base_traits_s2(const Constrained_delaunay_triangulation_base_traits_s2 & other) :
	m_traits(other.m_traits),
	m_apg(other.m_apg),
	m_proj(other.m_proj),
	m_significands(other.m_significands)
	{
		assert(m_significands > 2);
	}
	
	~Constrained_delaunay_triangulation_base_traits_s2() {}

	Side_of_oriented_circle_2 side_of_oriented_circle_2_object () const {
		return Side_of_oriented_circle_2(orientation_3_object());
	}
	
	Compare_distance_2 compare_distance_2_object () const {
		return baseTraits().compare_distance_3_object();
	}
	
	Orientation_2 orientation_2_object () const {
		return Orientation_2( orientation_3_object(), is_in_auxiliary_triangle_object());
	}
	
	Compare_x_2 compare_x_2_object() const {
		return Compare_x_2();
	}
	
	Compare_y_2 compare_y_2_object() const {
		return Compare_y_2();
	}
	
	Less_x_3 less_x_3_object() const {
		return baseTraits().less_x_3_object();
	}

	Less_y_3 less_y_3_object() const {
		return baseTraits().less_y_3_object();
	}

	Less_z_3 less_z_3_object() const {
		return baseTraits().less_z_3_object();
	}
	
	//stuff to support intersections
	Construct_segment_2 construct_segment_2_object() const {
		return construct_segment_3_object();
	}

	Construct_point_2 construct_point_2_object() const {
		return construct_point_3_object();
	}
	
	Less_x_2 less_x_2_object() const {
		return less_x_3_object();
	}

	Less_y_2 less_y_2_object() const {
		return less_y_3_object();
	}
	
	Project_on_sphere project_on_sphere_object() const {
		return Project_on_sphere(significands(), projector());
	}

	Project_on_sphere project_on_sphere_object(int _significands) const {
		return Project_on_sphere(_significands, projector());
	}
	
	Is_auxiliary_point is_auxiliary_point_object() const {
		return auxiliaryPointsGenerator().is_auxiliary_point_object();
	}
	
	Is_in_auxiliary_triangle is_in_auxiliary_triangle_object() const {
		return auxiliaryPointsGenerator().is_in_auxiliary_triangle_object();
	}
	
	Is_valid_point_on_sphere is_valid_point_on_sphere_object() const {
		return auxiliaryPointsGenerator().is_valid_point_on_sphere_object();
	}
	
	Generate_auxiliary_point generate_auxiliary_point_object() const {
		return auxiliaryPointsGenerator().generate_auxiliary_point_object();
	}
public:
	int significands() const { return m_significands; }
	const MyBaseTrait & baseTraits() const { return m_traits; }
	AuxiliaryPointsGenerator const & auxiliaryPointsGenerator() const { return m_apg; }
	const Projector & projector() const { return m_proj; }
protected:
	MyBaseTrait & baseTraits() { return m_traits; }
private:
	MyBaseTrait m_traits;
	AuxiliaryPointsGenerator m_apg;
	Projector m_proj;
	int m_significands;
};

#undef DEBUG_OUT

}//end namespace LIB_DTS2_NAMESPACE

#endif
