#include "gdal.h"
#include "parameters.h"
#include "merge_manager.h"


using namespace LxGeo::LxShapefilesMerge;


int main(int argc, char** argv) {
    
	clock_t t_begin = clock();
	GDALAllRegister();

	// Reads command-line parameters

	params = new Parameters(argc, argv);
	if (!params->initialized()) {
		delete params;
		return 1;
	}

	// Runs process

	MergeManager* mm = new MergeManager();
	mm->run();

	// Quits

	delete mm;
	delete params;

	clock_t t_end = clock();
	std::cout << "** Elapsed time : " << float(t_end - t_begin) / CLOCKS_PER_SEC << " s." << std::endl;

	return 0;
}
