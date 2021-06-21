#pragma once
#include "defs.h"



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

			void help();

			void parse(int argc, char* argv[]);

		public:
			bool printed_help;

			std::vector<std::string> paths;

			std::string output_basename;
			std::string output_shapefile;

			bool fix_srs_difference;
			double MAX_GROUPING_DISTANCE;
			double MAX_GROUPING_ANGLE_DEG;
			double e_distance_weight;
			double e_angle_weight;

			bool os_print;
			bool os_draw;
			bool os_check;
			double os_lambda_v;
			double os_lambda_c;
		};

		extern Parameters* params;
	}
}
