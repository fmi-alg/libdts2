#ifndef LIB_DTS2_CONSTRAINED_DELAUNAY_TRIANGULATION_BASE_TRAITS_S2_H
#define LIB_DTS2_CONSTRAINED_DELAUNAY_TRIANGULATION_BASE_TRAITS_S2_H

#include <libdts2/constants.h>
#include <libratss/ProjectS2.h>
#include <libratss/SphericalCoord.h>
#include <libratss/GeoCoord.h>

#include <CGAL/number_utils.h>

#include <CGAL/enum.h>


#include <assert.h>

namespace LIB_DTS2_NAMESPACE {

#ifndef NDEBUG
	#define DEBUG_OUT(__STR) {std::cout <<  __STR << std::endl;}
#else
	#define DEBUG_OUT(__) {}
#endif
#define LIB_DTS2_ORIGIN Point_3(0, 0, 0)

template<typename T_LINEAR_KERNEL>
class Constrained_delaunay_triangulation_base_traits_s2 {
public:
	using LinearKernel = T_LINEAR_KERNEL;
	using MyBaseTrait = T_LINEAR_KERNEL;
	using MyKernel = T_LINEAR_KERNEL;
	using MySelf = Constrained_delaunay_triangulation_base_traits_s2<MyBaseTrait>;
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
	using Point_3 = typename MyBaseTrait::Point_3;
	using Point_2 = Point_3;
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
		Orientation_2(const Orientation_3 & _ot3, const FT & _epsZ) : ot3(_ot3), m_epsZ(_epsZ) {}
		Orientation operator()(const Point_3 & p, const Point_3 & q, const Point_3 & r) const
		{
			assert(p.z() < 1 && q.z() < 1 && r.z() < 1);
			assert(p != q);
			
			Orientation oriented_side = ot3(p, q, LIB_DTS2_ORIGIN, r);
			
			if (oriented_side != Orientation::COLLINEAR && p.z() >= m_epsZ && q.z() >= m_epsZ) {
				//the point always has to be on the opposite site of the infinite vertex
				oriented_side = - ot3(p, q, LIB_DTS2_ORIGIN, Point(0, 0, 1));
			}
			return oriented_side;
		}
	private:
		Orientation_3 ot3;
		FT m_epsZ;
	};
	
	class Side_of_oriented_circle_2 {
	public:
		Side_of_oriented_circle_2(const Orientation_3 & _ot3) : ot3(_ot3) {}
		Orientation operator()(const Point_3 & p, const Point_3 & q, const Point_3 & r, const Point_3 & s) const
		{
			assert(p != q && p != r && q != r && r != s);
			auto oriented_side = ot3(p, q, r, s);
			return oriented_side;
		}
	private:
		Orientation_3 ot3;
	};
	
	class Compare_x_2 {
	public:
		Compare_x_2() {}
	public:
		Comparison_result  operator()(const Point_3 & p, const Point_3 & q) const {
			DEBUG_OUT("Compare_x_2 called")
			return m_cx3(p, q);
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
			DEBUG_OUT("Compare_y_2 called")
			//if this function gets called then we have to check why!
			assert(false);
			return m_cy3(p, q);
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
			if (!m_cop3(LIB_DTS2_ORIGIN, p, q, r)) {
				return false;
			}
			if (p == q || q == r) {
				return true;
			}
			//are coplanar, the direction of one point to another is defined by the vector with shorter euclidean distance
			//the following works since p, q, r have all the same length
			//here we check the angle >pq and >pr. >pq has to be smaller than >pr
			auto pq_sp = m_csp3(p, q);
			auto pr_sp = m_csp3(p, r);
			if (pq_sp > pr_sp) {
				return false;
			}
			//>pr has a larger angle than >pq BUT it can still be on the other other side of the circle
			//we have to make sure that r, q are on the same side of the line from p to -p
			
			//if the cross-product of pq and pr point the same direction,
			//then q,r are on the same side
			//They point in the same direction if pq_cpv = l*pr_cpv with l positive
			//They point in different direction if l negative
			auto pq_cpv = m_ccpv3(p, q);
			auto pr_cpv = m_ccpv3(p, r);
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
	public:
		//TODO:mark deprecated
		Point_3 operator()(const std::pair<double, double> & geoCoord) const {
			return this->operator()( LIB_RATSS_NAMESPACE::GeoCoord(geoCoord) );
		}
		Point_3 operator()(const LIB_RATSS_NAMESPACE::GeoCoord & geoCoord) const {
			mpq_class x, y, z;
			m_proj.projectFromGeo(geoCoord.lat, geoCoord.lon, x, y, z, significands());
			return Point_3(
				Conversion<FT>::moveFrom(x),
				Conversion<FT>::moveFrom(y),
				Conversion<FT>::moveFrom(z)
			);
		}
		Point_3 operator()(const LIB_RATSS_NAMESPACE::SphericalCoord & sphericalCoord) const {
			mpq_class x, y, z;
			m_proj.projectFromSpherical(sphericalCoord.theta, sphericalCoord.phi, x, y, z, significands());
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
		const Projector & projector() const {
			return m_proj;
		}
	private:
		int m_significands;
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
	Constrained_delaunay_triangulation_base_traits_s2() :
	m_epsZ(FT(1)-FT(std::numeric_limits<double>::epsilon()))
	{}
	///@param epsilon set the value of the z-coordinate above which no points should exist
	Constrained_delaunay_triangulation_base_traits_s2(const FT & _epsilon, int _significands) :
	m_epsZ(_epsilon),
	m_significands(_significands)
	{
		assert(m_significands > 2);
	}
	
	Constrained_delaunay_triangulation_base_traits_s2(const Constrained_delaunay_triangulation_base_traits_s2 & other) :
	m_traits(other.m_traits),
	m_epsZ(other.m_epsZ),
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
		return Orientation_2( orientation_3_object(), m_epsZ);
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
	
public:
	int significands() const { return m_significands; }
	const FT & epsZ() const { return m_epsZ; }
	const MyBaseTrait & baseTraits() const { return m_traits; }
	const Projector & projector() const { return m_proj; }
protected:
	MyBaseTrait & baseTraits() { return m_traits; } 
private:
	MyBaseTrait m_traits;
	FT m_epsZ;
	Projector m_proj;
	int m_significands;
};

#undef DEBUG_OUT
#undef LIB_DTS2_ORIGIN

}//end namespace LIB_DTS2_NAMESPACE

#endif
