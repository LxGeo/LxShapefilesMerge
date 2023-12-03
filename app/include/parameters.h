#pragma once
#include "defs.h"
#include <boost/filesystem.hpp>


namespace LxGeo
{
	namespace LxShapefilesMerge
	{
		class Parameters
		{
		public:
			Parameters(int argc, char* argv[]);

			~Parameters();

			bool initialized();

		protected:
			void init();

			void post_parse();

			void help();

			void parse(int argc, char* argv[]);

		public:
			bool printed_help;

			std::vector<std::string> paths;

			std::string output_basename;
			std::string output_shapefile;
			std::string temp_dir;

			bool fix_srs_difference;
			double MAX_GROUPING_DISTANCE;
			double MAX_GROUPING_ANGLE_DEG;
			double MIN_SEG_OVERLAP_RATIO;
			double e_distance_weight;
			double e_angle_weight;

			bool os_print;
			bool os_draw;
			bool os_check;
			double os_lambda_v;
			double os_lambda_c;
			std::vector<boost::filesystem::path> regularized_layers_path;

			std::vector<size_t> edge_diff_layer_weights;

		};

		extern Parameters* params;
	}
}
