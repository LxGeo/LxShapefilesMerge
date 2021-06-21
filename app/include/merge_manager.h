#pragma once
#include "defs.h"

namespace LxGeo
{
	namespace LxShapefilesMerge
	{
		/**
		*  A MergeManager class to manage running required steps to generate final merge shapefile.
		*/
		class MergeManager
		{
		public:
			MergeManager();

			~MergeManager();

			/**
			*  A method used to run all steps of merging.
			*/
			virtual void run();

			/**
			*  A method to check requirements before running merging steps.
			* Example: -Checking shapefiles before loading, -Checking shapefile layers share the same srs, ...
			* @return an int indicating error severity
			*/
			int pre_check();

		public:
			std::vector<std::string>& all_input_shapefiles;
			std::string& output_shapefile;
			bool apply_srs_transform = false;
		};
	}
}