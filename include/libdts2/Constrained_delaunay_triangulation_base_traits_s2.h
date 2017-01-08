#ifndef LIB_DTS2_CONSTRAINED_DELAUNAY_TRIANGULATION_BASE_TRAITS_S2_H
#define LIB_DTS2_CONSTRAINED_DELAUNAY_TRIANGULATION_BASE_TRAITS_S2_H

#include <libdts2/constants.h>
#include <libratss/ProjectS2.h>
#include <libratss/SphericalCoord.h>
#include <libratss/GeoCoord.h>

#include <CGAL/enum.h>


#include <assert.h>

namespace LIB_DTS2_NAMESPACE {

#define DEBUG_OUT(__STR) {std::cout <<  __STR << std::endl;}
#define LIB_DTS2_ORIGIN Point_3(0, 0, 0)

template<typename T_KERNEL>
class Constrained_delaunay_triangulation_base_traits_s2 {
public:
	using MyBaseTrait = T_KERNEL;
	using MyKernel = T_KERNEL;
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
	using Intersect_3 = typename MyBaseTrait::Intersect_3;
	using Plane_3 = typename MyBaseTrait::Plane_3;
	using Line_3 = typename MyBaseTrait::Line_3;
	using Ray_3 = typename MyBaseTrait::Ray_3;
	using Vector_3 = typename MyBaseTrait::Vector_3;
	
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
			
			Orientation oriented_side;
			
			if (p.z() >= m_epsZ && q.z() >= m_epsZ) {
				//the point always has to be on the opposite site of the infinite vertex
				oriented_side = - ot3(p, q, LIB_DTS2_ORIGIN, Point(0, 0, 1));
			}
			else {
				oriented_side = ot3(p, q, LIB_DTS2_ORIGIN, r);
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
public:
	
	///sub classes of this trait should sub-class this class and add an appropriate
	///Point_3 operator()(const Point_3 & p) const; function to map points to the sphere 
	class Project_on_sphere {
	private:
		template<typename T>
		using Conversion = LIB_RATSS_NAMESPACE::Conversion<T>;
	public:
		Project_on_sphere(int _precision, const Projector & p) : m_precision(_precision), m_proj(p) {}
	public:
		//TODO:mark deprecated
		Point_3 operator()(const std::pair<double, double> & geoCoord) const {
			return this->operator()( LIB_RATSS_NAMESPACE::GeoCoord(geoCoord) );
		}
		Point_3 operator()(const LIB_RATSS_NAMESPACE::GeoCoord & geoCoord) const {
			mpq_class x, y, z;
			m_proj.projectFromGeo(geoCoord.lat, geoCoord.lon, x, y, z, precision());
			return Point_3(
				Conversion<FT>::moveFrom(x),
				Conversion<FT>::moveFrom(y),
				Conversion<FT>::moveFrom(z)
			);
		}
		Point_3 operator()(const LIB_RATSS_NAMESPACE::SphericalCoord & sphericalCoord) const {
			mpq_class x, y, z;
			m_proj.projectFromSpherical(sphericalCoord.theta, sphericalCoord.phi, x, y, z, precision());
			return Point_3(
				Conversion<FT>::moveFrom(x),
				Conversion<FT>::moveFrom(y),
				Conversion<FT>::moveFrom(z)
			);
		}
	public:
		int precision() const {
			return m_precision;
		}
		const Projector & projector() const {
			return m_proj;
		}
	private:
		int m_precision;
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
protected:
	static inline void print(const Point_3 & p3) {
		std::cout << '(' << doubleValue(p3.x()) << ','
		<< doubleValue(p3.y()) << ','
		<< doubleValue(p3.z()) << ')' << '\n';
	}
	static inline void print(std::ostream & out, const Point_3 & p) {
		using Conversion = LIB_RATSS_NAMESPACE::Conversion<FT>;
		out << '(' << Conversion::toMpq(p.x()) << ", " << Conversion::toMpq(p.y()) << ", " << Conversion::toMpq(p.z()) << ')';
	}
public:
	const FT & epsilon() const;
public:
	static inline double doubleValue(const FT & v) {
		return LIB_RATSS_NAMESPACE::Conversion<FT>::toMpreal(v, 53).toDouble();
	}
public: //object functions
	///This does not correctly initialize this trait!
	Constrained_delaunay_triangulation_base_traits_s2() :
	m_epsZ(FT(std::numeric_limits<double>::epsilon()))
	{}
	///@param epsilon set the value of the z-coordinate above which no points should exist
	Constrained_delaunay_triangulation_base_traits_s2(const FT & _epsilon, int _precision) :
	m_epsZ(_epsilon),
	m_precision(_precision)
	{}
	
	Constrained_delaunay_triangulation_base_traits_s2(const Constrained_delaunay_triangulation_base_traits_s2 & other) :
	m_traits(other.m_traits),
	m_epsZ(other.m_epsZ),
	m_proj(other.m_proj),
	m_precision(other.m_precision)
	{}
	
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
	
	Less_x_2 less_x_2_object() const {
		return less_x_3_object();
	}

	Less_y_2 less_y_2_object() const {
		return less_y_3_object();
	}
	
	Project_on_sphere project_on_sphere_object() const {
		return Project_on_sphere(precision(), projector());
	}
	
public:
	int precision() const { return m_precision; }
	const FT & epsZ() const { return m_epsZ; }
	const MyBaseTrait & baseTraits() const { return m_traits; }
	const Projector & projector() const { return m_proj; }
protected:
	MyBaseTrait & baseTraits() { return m_traits; } 
private:
	MyBaseTrait m_traits;
	FT m_epsZ;
	Projector m_proj;
	int m_precision;
};

#undef DEBUG_OUT
#undef LIB_DTS2_ORIGIN

}//end namespace LIB_DTS2_NAMESPACE

#endif