#pragma once
#pragma warning( disable : 4251 )

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

namespace LxGeo
{
	namespace GeometryFactoryShared
	{
		namespace bg = boost::geometry;
		namespace bgi = boost::geometry::index;
		typedef bg::model::point<double, 2, bg::cs::cartesian> Boost_Point_2;
		typedef bg::model::box<Boost_Point_2> Boost_Box_2;
		typedef std::pair<Boost_Box_2, size_t> Boost_Value_2;
		typedef bgi::rtree<Boost_Value_2, bgi::quadratic<16> > Boost_RTree_2;
		typedef bg::model::segment<Boost_Point_2> Boost_Segment_2;
	}
}