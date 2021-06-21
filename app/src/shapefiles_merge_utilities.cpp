#include "shapefiles_merge_utilities.h"
#include "parameters.h"
#include "defs.h"
#include "defs_cgal.h"
#include "defs_boost.h"
#include <gdal.h>
#include <gdal_priv.h>
#include <ogrsf_frmts.h>
#include "segment_graph.h"

namespace LxGeo
{
	namespace LxShapefilesMerge
	{
		namespace shapefiles_merge_utils
		{

			using namespace GeometryFactoryShared;

			int check_shapefiles_validity(std::vector<std::string>& all_paths, bool& apply_srs_transform)
			{
				BOOST_LOG_TRIVIAL(info) << "Checking shapefiles validity!";
				int validation_severity = NO_ERROR;
				GDALDataset* dataset = NULL;
				std::set<OGRSpatialReference> srs_set;

				try {
					for (size_t i = 0; i < all_paths.size(); ++i) {

						const char* c_dataset_path = all_paths[i].c_str();
						dataset = (GDALDataset*)GDALOpenEx(c_dataset_path, GDAL_OF_VECTOR, NULL, NULL, NULL);
						if (dataset == NULL) {
							BOOST_LOG_TRIVIAL(warning) << "Input shapefile path: " << c_dataset_path << " is missing";
							validation_severity = (validation_severity > WRONG_SHAPEFILE_PATH) ? validation_severity : WRONG_SHAPEFILE_PATH;
						}
						else
						{
							OGRLayer* layer = dataset->GetLayer(0);
							OGRSpatialReference* source_srs = layer->GetSpatialRef();
							srs_set.insert(*source_srs);
							if (srs_set.size() > 1)
							{
								BOOST_LOG_TRIVIAL(warning) << "One of the input shapefiles do not share the same SpatialRefSystem!";
								validation_severity = (validation_severity > SPATIAL_REF_CONFLICT) ? validation_severity : SPATIAL_REF_CONFLICT;
								BOOST_LOG_TRIVIAL(info) << "Applying srs transformation based on the first shapefile srs!";
								if (params->fix_srs_difference) apply_srs_transform = true;
							}
						}
					}
				}
				catch (std::exception& e) {
					BOOST_LOG_TRIVIAL(debug) << e.what();
					BOOST_LOG_TRIVIAL(fatal) << "Fatal error! For a quick resolution keep a copy of input data!";
					validation_severity = UNKNOWN_ERROR;
				}

				return validation_severity;

			}

			void load_segments_data(std::vector<std::string>& all_paths,
				std::vector<Inexact_Segment_2>& all_segments,
				std::vector<short int>& segment_LID,
				std::vector<short int>& segment_PID,
				std::vector<short int>& segment_ORDinP,
				const bool apply_srs_transform)
			{
				if (all_paths.empty()) return;

				GDALDataset* dataset_ref = NULL;
				GDALDataset* dataset = NULL;

				try {

					// Gets a first spatial reference
					dataset_ref = (GDALDataset*)GDALOpenEx(all_paths[0].c_str(), GDAL_OF_VECTOR, NULL, NULL, NULL);
					if (dataset_ref == NULL) {
						throw std::ios_base::failure("Error : unable to load shapefile.");
					}
					OGRLayer* layer_ref = dataset_ref->GetLayer(0);
					OGRSpatialReference* target_ref = layer_ref->GetSpatialRef();

					// reserve segments vector size = ( #of_shapefile * #of_polygons * mean_of_segments_per_polygon)
					int mean_of_segments_per_polygon = 6;
					size_t estimated_segments_count = all_paths.size() * layer_ref->GetFeatureCount() * mean_of_segments_per_polygon;
					BOOST_LOG_TRIVIAL(info) << "estimated_segments_count :" << estimated_segments_count;
					all_segments.reserve(estimated_segments_count);
					segment_LID.reserve(estimated_segments_count);
					segment_PID.reserve(estimated_segments_count);
					segment_ORDinP.reserve(estimated_segments_count);


					// iterate through shapefiles
					for (size_t i = 0; i < all_paths.size(); ++i) {
						GDALDataset* dataset = NULL;

						// Step 1.
						// Opens the i-th shapefile
						// Creates a coordinate transformation towards the reference coordinate system 
						if (i == 0) {
							dataset = dataset_ref;
						}
						else {
							dataset = (GDALDataset*)GDALOpenEx(all_paths[i].c_str(), GDAL_OF_VECTOR, NULL, NULL, NULL);
						}

						if (dataset == NULL) {
							throw std::ios_base::failure("Error : unable to load shapefile.");
						}

						OGRLayer* layer = dataset->GetLayer(0);
						OGRSpatialReference* source_ref = layer->GetSpatialRef();
						OGRCoordinateTransformation* CT = OGRCreateCoordinateTransformation(source_ref, target_ref);

						// Step 2.
						// Reads contents of the current shapefile

						size_t n = size_t(layer->GetFeatureCount());

						for (size_t j = 0; j < n; ++j) {
							OGRFeature* feat = layer->GetNextFeature();
							if (feat == NULL) continue;

							OGRGeometry* geom = feat->GetGeometryRef();

							// Assumes the shapefile only contains OGRPolygons
							if (OGRPolygon* P = dynamic_cast<OGRPolygon*>(geom)) {

								// Reads rings
								std::list<OGRLinearRing*> ogr_rings;
								ogr_rings.push_back(P->getExteriorRing());
								for (int k = 0; k < P->getNumInteriorRings(); ++k) ogr_rings.push_back(P->getInteriorRing(k));

								for (OGRLinearRing* ring : ogr_rings) {
									
									// applying srs transformation if requested
									std::vector<Inexact_Point_2> R;
									int ring_size = ring->getNumPoints();
									R.reserve(ring_size);
									for (int u = 0; u < ring_size - 1; ++u) {
										OGRPoint pt;
										ring->getPoint(u, &pt);
										double x = pt.getX(), y = pt.getY();
										if (apply_srs_transform) CT->Transform(1, &x, &y);
										R.push_back(Inexact_Point_2(x, y));
									}

									// Creating segments
									for (int l = 0; l < ring_size-2; ++l) {
										Inexact_Point_2 p1 = R[l];
										Inexact_Point_2 p2 = R[l+1];
										// adding to segments
										all_segments.push_back(Inexact_Segment_2(p1, p2));
										segment_LID.push_back(i);
										segment_PID.push_back(j);
										// if two segments share the same ORDinP & PID than at least an outer ring exists.
										segment_ORDinP.push_back(l);
									}

								}

							}
						}

						// Step 3.
						// Reading is over, closes the dataset
						if (i != 0) GDALClose(dataset);
					}

				}
				catch (std::exception& e) {
					BOOST_LOG_TRIVIAL(debug) << e.what();
					BOOST_LOG_TRIVIAL(fatal) << "Fatal error! For a quick resolution keep a copy of input data!";

					all_segments.clear();
					segment_LID.clear();
					segment_PID.clear();
					segment_ORDinP.clear();
				}

				if (dataset_ref != NULL) GDALClose(dataset_ref);
				if (dataset != NULL) GDALClose(dataset);

			}

			void regularize_segments(std::vector<Inexact_Segment_2>& all_segments,
				std::vector<short int>& segment_LID,
				std::vector<short int>& segment_PID,
				std::vector<short int>& segment_ORDinP) {

				BOOST_LOG_TRIVIAL(info) << "Computing segments angles!";
				std::vector<double> segment_angle;
				segment_angle.reserve(all_segments.size());
				for (auto& c_segment : all_segments) {
					Inexact_Point_2 p1 = c_segment.source();
					Inexact_Point_2 p2 = c_segment.target();
					double atan_val = std::atan2(p2.y() - p1.y(), p2.x() - p1.x());
					segment_angle.push_back(fmod(atan_val, M_PI));
				};

				BOOST_LOG_TRIVIAL(info) << "Computing segments tree!";
				Boost_RTree_2 segments_tree;
				clock_t t_begin_segments_tree = clock();
				make_rtree(all_segments, segments_tree);
				clock_t t_end_segments_tree = clock();
				BOOST_LOG_TRIVIAL(debug) << "Elapsed time for segments tree creation: " << double(t_end_segments_tree - t_begin_segments_tree) / CLOCKS_PER_SEC << " s.";
				
				SegmentGraph* SG = new SegmentGraph(all_segments,
					segment_LID, segment_PID, segment_angle, segments_tree);


			}



			void make_rtree(const std::vector< Inexact_Segment_2 >& all_segments, Boost_RTree_2& RT)
			{
				for (size_t i = 0; i < all_segments.size(); ++i) {
					const Inexact_Segment_2 c_segment = all_segments[i];

					// Computes bbox of each segment

					Bbox_2 B = CGAL::bbox_2(c_segment.source(), c_segment.target());
					double xmin = B.xmin(), xmax = B.xmax(), ymin = B.ymin(), ymax = B.ymax();

					Boost_Box_2 box(Boost_Point_2(xmin, ymin), Boost_Point_2(xmax, ymax));
					RT.insert(Boost_Value_2(box, i));
				}
			}

			
			
			//void read_shapefiles(const std::vector<std::string>& all_paths, std::vector<std::vector<std::vector<std::vector<Inexact_Point_2> > > >& all_polygons);

			//void recenter(std::vector<std::vector<std::vector<std::vector<Inexact_Point_2> > > >& all_polygons, Inexact_Vector_2& shift, Bbox_2& bbox);

			//void make_segments(const std::vector<std::vector<std::vector<std::vector<Inexact_Point_2> > > >& all_polygons, std::vector<Inexact_Segment_2>& all_segments);

			//void write_shapefile(const std::string& input_filename, const std::string& output_filename, const std::vector<Polygon_with_attributes*>& polygons);
		}
	}
}