#pragma once
#include "defs.h"
#include <boost/filesystem.hpp>
#include "defs_cgal.h"
#include "defs_boost.h"
#include <gdal.h>
#include <gdal_priv.h>

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

			std::string overlay_union_layers(std::vector<boost::filesystem::path> &regularized_layers_path);

			double segments_overlap_ratios(const Inexact_Segment_2& segment1, const Inexact_Segment_2& segment2);

			void clean_invalid(OGRLayer* c_layer);

			void simplify_ring(std::vector<Inexact_Point_2>& R, std::vector<Inexact_Point_2>& simplified_R);

			void read_shapefiles(const std::vector<std::string>& all_paths, std::vector<std::vector<std::vector<std::vector<Inexact_Point_2> > > >& all_polygons);

			void extract_edges_from_polygons(std::vector<std::vector<std::vector<Inexact_Point_2>>>& c_layer_regularized_polygons,
				std::vector<Boost_LineString_2>& c_layer_regularized_edges);

			void make_rtree_polygons(const std::vector<std::vector<std::vector<Inexact_Point_2> > >& polygons, Boost_RTree_2& RT);

			void make_rtree_linestrings(const std::vector<Boost_LineString_2>& linestrings, Boost_RTree_2& RT);

			void read_single_shapefile(std::string shapefile_path, OGRSpatialReference* target_ref, std::vector<std::vector<std::vector<Inexact_Point_2> > >& polygons);

			void read_single_shapefile(GDALDataset* shapefile_dataset, OGRSpatialReference* target_ref, std::vector<std::vector<std::vector<Inexact_Point_2> > >& polygons);

		}
	}
}