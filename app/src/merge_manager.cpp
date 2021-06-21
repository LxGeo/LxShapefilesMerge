#include "merge_manager.h"
#include "defs_cgal.h"
#include "parameters.h"
#include "shapefiles_merge_utilities.h"


namespace LxGeo
{
	namespace LxShapefilesMerge
	{
		MergeManager::MergeManager()
		{
			all_input_shapefiles = params->paths;
			output_shapefile = params->output_shapefile;

		}


		MergeManager::~MergeManager()
		{
		}

		int MergeManager::pre_check() {

			using namespace shapefiles_merge_utils;

			return check_shapefiles_validity(all_input_shapefiles, apply_srs_transform);

		}

		void MergeManager::run()
		{
			using namespace shapefiles_merge_utils;
			using namespace GeometryFactoryShared;

			// Checking before ruuning pipeline
			int validation_severity = pre_check();
			if (validation_severity >= WRONG_SHAPEFILE_PATH) return;
			else if (validation_severity == SPATIAL_REF_CONFLICT && !params->fix_srs_difference) return;
			else BOOST_LOG_TRIVIAL(info) << "Pre check passed succefully!";

			// Loading polygons as segments (to run regularization)

			std::vector<Inexact_Segment_2> all_segments;
			// segment_LID contains segments respective layer ID
			std::vector<short int> segment_LID;
			// segment_PID contains segments respective polygon ID
			std::vector<short int> segment_PID;
			// segment_ORDinP contains segments respective segment order in a polygon
			std::vector<short int> segment_ORDinP;

			load_segments_data(all_input_shapefiles,
				all_segments,
				segment_LID,
				segment_PID,
				segment_LID,
				apply_srs_transform);

			// Regularizing Segments
			regularize_segments(all_segments,
				segment_LID,
				segment_PID,
				segment_LID);


		}
	}
}