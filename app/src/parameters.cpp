#include "parameters.h"

namespace LxGeo
{
	namespace LxShapefilesMerge
	{
		Parameters::Parameters(int argc, char* argv[])
		{
			init();
			parse(argc, argv);
			
		}


		Parameters::~Parameters()
		{
		}


		bool Parameters::initialized()
		{
			bool is_initialized = !paths.empty();
			if (!is_initialized) help();

			return is_initialized;
		}


		void Parameters::init()
		{
			printed_help = false;

			paths.clear();

			output_basename = "result";
			output_shapefile = "result.shp";
			temp_dir = "temp_dir";


			MAX_GROUPING_DISTANCE = 2;
			MAX_GROUPING_ANGLE_DEG = 3;
			e_distance_weight=0.5;
			e_angle_weight=0.5;
			fix_srs_difference = os_print = os_draw = os_check = false;
			os_lambda_v = 0.001;
			os_lambda_c = 1;
		}


		void Parameters::help()
		{
			if (printed_help) return;

			std::cout << "lxShapefileFusion.exe [args]" << std::endl
				<< std::endl
				<< "where [args] are : " << std::endl
				<< "  [-h | --help] -> print this help message" << std::endl
				<< std::endl
				<< "** Basic parameters : " << std::endl
				<< std::endl
				<< "  [-i] [n] [path_1] [path_2] ... [path_n] -> provide paths of input shapefiles" << std::endl
				<< "  [-o] [basename] -> specify basename of output file" << std::endl
				<< "  [--fix_srs_difference] to apply transformation of SpatialRefrenceSystem if any is different." << std::endl
				<< "  [-mg_distance] [MAX_GROUPING_DISTANCE] -> specify maximum grouping distance for segments (meters). Default 2." << std::endl
				<< "  [-mg_angle] [MAX_GROUPING_ANGLE_DEG] -> specify maximum grouping angle for segments (degrees). Default 3." << std::endl
				<< std::endl
				<< "Version compiled on : " << __DATE__ << std::endl;

			printed_help = true;
		}


		void Parameters::parse(int argc, char* argv[])
		{
			if (argc == 1) {
				help();
				return;
			}

			std::list<std::string> unknown_args;

			size_t r = 1;
			while (r < argc) {
				std::string arg = argv[r];
				if (arg == "-h" || arg == "--help") {
					help();
					return;
				}
				else if (arg == "-i" && r + 1 < argc) {
					size_t n = atoi(argv[r + 1]);
					if (r + 1 + n < argc) {
						for (size_t i = 0; i < n; ++i) paths.push_back(argv[r + 2 + i]);
					}
					r += (n + 2);
				}
				else if ((arg == "-o" || arg == "--output") && r + 1 < argc) {
					std::string f = argv[r + 1];
					std::string extension = (f.size() > 4 ? f.substr(f.size() - 4, 4) : "");
					if (extension == ".shp" || extension == ".SHP") {
						output_shapefile = f;
						output_basename = f.substr(0, f.size() - 4);
					}
					else {
						std::cout << "Warning : invalid output filename. Writing result in result.shp" << std::endl;
					}
					r += 2;

				}
				else if (arg == "-mg_distance") {
					MAX_GROUPING_DISTANCE = atof(argv[r + 1]);
					r += 1;
				}
				else if (arg == "-mg_angle") {
					MAX_GROUPING_ANGLE_DEG = atof(argv[r + 1]);
					r += 1;
				}
				else if (arg == "--fix_srs_difference") {
					fix_srs_difference = true;
					r += 1;
				}
				else if (arg == "--os-print") {
					os_print = true;
					r += 1;
				}
				else if (arg == "--os-draw") {
					os_draw = true;
					r += 1;
				}
				else if (arg == "--os-check") {
					os_check = true;
					r += 1;
				}
				else if (arg == "--osl" && r + 2 < argc) {
					os_lambda_v = atof(argv[r + 1]);
					os_lambda_c = atof(argv[r + 2]);
					r += 3;
				}
				else {
					unknown_args.push_back(arg);
					r += 1;
				}
			}

			if (!unknown_args.empty()) {
				std::cout << "There were unknown arguments in command line call :" << std::endl << '\t';
				for (const std::string& arg : unknown_args) std::cout << arg << " ";
				std::cout << std::endl;
				help();
			}
		}

		Parameters* params = nullptr;
	}
}