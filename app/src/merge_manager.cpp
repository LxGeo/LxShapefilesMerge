#include "merge_manager.h"
#include "defs_cgal.h"
#include "parameters.h"
#include "shapefiles_merge_utilities.h"
#include "overlay_labeler.h"


namespace LxGeo
{
	namespace LxShapefilesMerge
	{
		MergeManager::MergeManager() :all_input_shapefiles(params->paths), output_shapefile(params->output_shapefile)
		{
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
			if (validation_severity >= S_WRONG_SHAPEFILE_PATH) return;
			else if ((validation_severity == S_SPATIAL_REF_CONFLICT) && !params->fix_srs_difference) return;
			else if (validation_severity == S_DIRECTORY_CREATION_ERROR)
			{
				BOOST_LOG_TRIVIAL(fatal) << "Cannot start merging steps! Ensure access to output directory and temporary directory!";
				return;
			}
			else BOOST_LOG_TRIVIAL(info) << "Pre check passed succefully!";

			// Loading polygons as segments (to run regularization)

			std::vector<Segment_2> all_segments;
			// segment_LID contains segments respective layer ID
			std::vector<short int> segment_LID;
			// segment_PID contains segments respective polygon ID
			std::vector<short int> segment_PID;
			// segment_ORDinP contains segments respective segment order in a polygon
			std::vector<short int> segment_ORDinP;
			// segment_resp_ring_size
			std::vector<short int> segment_RRSize;

			load_segments_data(all_input_shapefiles,
				all_segments,
				segment_LID,
				segment_PID,
				segment_ORDinP,
				segment_RRSize,
				apply_srs_transform);

			// file_path constructions
			for (size_t layer_id = 0; layer_id < params->paths.size(); ++layer_id) {
				boost::filesystem::path c_out_basename(fmt::format("l_{}.shp", layer_id));
				boost::filesystem::path outdir(params->temp_dir);
				boost::filesystem::path c_full_path = outdir / c_out_basename;
				params->regularized_layers_path.push_back(c_full_path);
			}

			// Regularizing Segments
			regularize_segments(all_segments,
				segment_LID,
				segment_PID,
				segment_ORDinP,
				segment_RRSize);


			std::string total_overlay_layer_path;
			//total_overlay_layer_path = overlay_union_layers(params->regularized_layers_path);

			all_segments.clear();
			segment_LID.clear();
			segment_PID.clear();
			segment_ORDinP.clear();
			
			return;

			// total overlay union path
			std::string temp_overlay_layer_path = params->temp_dir + "\\grass_fusion\\U_0_1_2_3.shp";
			//  facet labeling optimizer

			std::vector<std::string> regularized_layers_path_string;
			for (auto reg_path : params->regularized_layers_path) regularized_layers_path_string.push_back(reg_path.string());
			OverlayLabeler Ol = OverlayLabeler(temp_overlay_layer_path,
				regularized_layers_path_string);

			Ol.construct_graph();

			Ol.log_subgraphs();

		}
	}
}