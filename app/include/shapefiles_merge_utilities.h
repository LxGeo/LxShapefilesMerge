#pragma once
#include "defs.h"
#include "defs_cgal.h"

namespace LxGeo
{
	namespace LxShapefilesMerge
	{
		namespace shapefiles_merge_utils
		{
			using namespace GeometryFactoryShared;
			

			/**
			*  A method to check requirements before running merging steps.
			*  @param all_paths: all shapefiles paths
			*  @return an int indicating error severity where (0: no errors; 1: at least one shapefile have a different SpatialRef; 2: at least one shapefile is missing
				3: Unknown error).
			*/
			int check_shapefiles_validity(std::vector<std::string>& all_paths, bool& apply_srs_transform);

			void load_segments_data(std::vector<std::string>& all_paths,
				std::vector<Inexact_Segment_2>& all_segments,
				std::vector<short int>& segment_LID,
				std::vector<short int>& segment_PID,
				std::vector<short int>& segment_ORDinP,
				const bool apply_srs_transform);

			void regularize_segments(std::vector<Inexact_Segment_2>& all_segments,
				std::vector<short int>& segment_LID,
				std::vector<short int>& segment_PID,
				std::vector<short int>& segment_ORDinP);

			
			//void read_shapefiles(const std::vector<std::string>& all_paths, std::vector<std::vector<std::vector<std::vector<Inexact_Point_2> > > >& all_polygons);

			//void recenter(std::vector<std::vector<std::vector<std::vector<Inexact_Point_2> > > >& all_polygons, Inexact_Vector_2& shift, Bbox_2& bbox);

			//void make_segments(const std::vector<std::vector<std::vector<std::vector<Inexact_Point_2> > > >& all_polygons, std::vector<Inexact_Segment_2>& all_segments);

			//void write_shapefile(const std::string& input_filename, const std::string& output_filename, const std::vector<Polygon_with_attributes*>& polygons);
		}
	}
}