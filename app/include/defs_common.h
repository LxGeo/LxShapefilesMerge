#pragma once

#include "defs_boost.h"
#include "defs_cgal.h"

namespace LxGeo
{
	namespace GeometryFactoryShared
	{

		Inexact_Point_2 transform_B2C_Point(const Boost_Point_2& boost_point) {
			return Inexact_Point_2(bg::get<0>(boost_point), bg::get<1>(boost_point));
		};

		Boost_Point_2 transform_C2B_Point(Inexact_Point_2& cgal_point) {
			return Boost_Point_2(cgal_point.x(), cgal_point.y());
		};

		void container_transform_B2C_Points(std::vector<Boost_Point_2>& input_container, std::vector<Inexact_Point_2>& output_container) {
			assert(output_container.empty());
			output_container.reserve(input_container.size());

			for (Boost_Point_2 c_point : input_container) {
				Inexact_Point_2 cgal_point = transform_B2C_Point(c_point);
				output_container.push_back(cgal_point);
			}
		};
	}
}