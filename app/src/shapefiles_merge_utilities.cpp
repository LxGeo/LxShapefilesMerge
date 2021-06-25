#include "shapefiles_merge_utilities.h"
#include "parameters.h"
#include "defs.h"
#include "defs_cgal.h"
#include "defs_boost.h"
#include <gdal.h>
#include <gdal_priv.h>
#include <ogrsf_frmts.h>
#include "segment_graph.h"
#include <boost/filesystem.hpp>
#include <fmt/core.h>

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
				int validation_severity = S_NO_ERROR;
				GDALDataset* dataset = NULL;
				std::set<int> srs_set;

				try {
					for (size_t i = 0; i < all_paths.size(); ++i) {

						const char* c_dataset_path = all_paths[i].c_str();
						dataset = (GDALDataset*)GDALOpenEx(c_dataset_path, GDAL_OF_VECTOR, NULL, NULL, NULL);
						if (dataset == NULL) {
							BOOST_LOG_TRIVIAL(warning) << "Input shapefile path: " << c_dataset_path << " is missing";
							validation_severity = (validation_severity > S_WRONG_SHAPEFILE_PATH) ? validation_severity : S_WRONG_SHAPEFILE_PATH;
						}
						else
						{
							OGRLayer* layer = dataset->GetLayer(0);
							OGRSpatialReference* source_srs = layer->GetSpatialRef();
							int utm_zone_number = source_srs->GetUTMZone();
							srs_set.insert(utm_zone_number);
														
						}
					}
					if (srs_set.size() > 1)
					{
						BOOST_LOG_TRIVIAL(warning) << "One of the input shapefiles do not share the same SpatialRefSystem!";
						validation_severity = (validation_severity > S_SPATIAL_REF_CONFLICT) ? validation_severity : S_SPATIAL_REF_CONFLICT;
						
						if (params->fix_srs_difference)
						{
							apply_srs_transform = true;
							BOOST_LOG_TRIVIAL(info) << "Applying srs transformation based on the first shapefile srs!";
						}
						else {
							BOOST_LOG_TRIVIAL(warning) << "Try setting --fix_srs_difference flag to apply srs tranformation!";
						}
					}

					//output dirs creation
					boost::filesystem::path output_path(params->output_shapefile);
					boost::filesystem::path output_parent_dirname = output_path.parent_path();
					boost::filesystem::path output_temp_path = output_parent_dirname / params->temp_dir;
					params->temp_dir = output_temp_path.string();
					if (boost::filesystem::exists(output_parent_dirname) || boost::filesystem::create_directory(output_parent_dirname))
					{
						BOOST_LOG_TRIVIAL(info) << fmt::format("Directory Created: {}", output_parent_dirname.string());
					}
					else validation_severity = S_DIRECTORY_CREATION_ERROR;
					if (boost::filesystem::exists(output_temp_path) ||  boost::filesystem::create_directory(output_temp_path))
					{
						BOOST_LOG_TRIVIAL(info) << fmt::format("Directory Created: {}", output_temp_path.string());
					}
					else validation_severity = S_DIRECTORY_CREATION_ERROR;

				}
				catch (std::exception& e) {
					BOOST_LOG_TRIVIAL(debug) << e.what();
					BOOST_LOG_TRIVIAL(fatal) << "Fatal error! For a quick resolution keep a copy of input data!";
					validation_severity = S_UNKNOWN_ERROR;
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
									for (int u = 0; u < ring_size; ++u) {
										OGRPoint pt;
										ring->getPoint(u, &pt);
										double x = pt.getX(), y = pt.getY();
										if (apply_srs_transform) CT->Transform(1, &x, &y);
										R.push_back(Inexact_Point_2(x, y));
									}

									// Creating segments
									size_t added_segments = 0;
									for (int l = 0; l < ring_size-2; ) {
										Inexact_Point_2 p1 = R[l];
										Inexact_Point_2 p2 = R[l+1];
										Inexact_Point_2 p3 = R[l + 2];

										std::vector<Inexact_Segment_2> segments_to_add(2);
										//calc_angle
										double numerator = p2.y() * (p1.x() - p3.x()) + p1.y() * (p3.x() - p2.x()) + p3.y() * (p2.x() - p1.x());
										double denominator = (p2.x() - p1.x()) * (p1.x() - p3.x()) + (p2.y() - p1.y()) * (p1.y() - p3.y());
										double ratio = numerator / denominator;
										double points_angle = std::abs(std::fmod(std::atan(ratio), M_PI));
										bool pts_collinear = (points_angle < ( M_PI / 180)); //CGAL::collinear(p1, p2, p3);
										if (pts_collinear) {
											segments_to_add.push_back(Inexact_Segment_2(p1, p3));
											l += 2;
										}
										else {
											segments_to_add.push_back(Inexact_Segment_2(p1, p2));
											l += 1;
											if (l== ring_size-2) segments_to_add.push_back(Inexact_Segment_2(p2, p3));
										}

										for (const auto& segment_to_add : segments_to_add) {
											if (segment_to_add.squared_length() == 0)
											{
												BOOST_LOG_TRIVIAL(debug) << "Ignoring segment for 0 length!";
												continue;
											}
											// adding to segments
											all_segments.push_back(segment_to_add);
											segment_LID.push_back(i);
											segment_PID.push_back(j);
											// if two segments share the same ORDinP & PID than at least an outer ring exists.
											segment_ORDinP.push_back(added_segments);
											added_segments++;
										}
										
									}

								}

							}
						}

						// Step 3.
						// Reading is over, closes the dataset
						if (i != 0) GDALClose(dataset);

						all_segments.shrink_to_fit();
						segment_LID.shrink_to_fit();
						segment_PID.shrink_to_fit();
						segment_ORDinP.shrink_to_fit();
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


			void make_rtree(const std::vector< Inexact_Segment_2 >& all_segments, Boost_RTree_2& RT)
			{
				for (size_t i = 0; i < all_segments.size(); ++i) {
					const Inexact_Segment_2 c_segment = all_segments[i];
										
					double xmin = std::min<double>(c_segment.source().x(), c_segment.target().x());
					double xmax = std::max<double>(c_segment.source().x(), c_segment.target().x());
					double ymin = std::min<double>(c_segment.source().y(), c_segment.target().y());
					double ymax = std::max<double>(c_segment.source().y(), c_segment.target().y());

					Boost_Box_2 box(Boost_Point_2(xmin, ymin), Boost_Point_2(xmax, ymax));
					RT.insert(Boost_Value_2(box, i));
				}
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
				BOOST_LOG_TRIVIAL(debug) << "Elapsed time for segments tree creation: " << float(t_end_segments_tree - t_begin_segments_tree) / CLOCKS_PER_SEC << " s.";
				
				SegmentGraph* SG = new SegmentGraph(all_segments,
					segment_LID, segment_PID, segment_ORDinP, segment_angle, segments_tree);

				//SG->write_grouped_segments_shapefile(params->output_shapefile);
				SG->fill_graph();
				SG->cluster_segments();
				SG->fuse_segments();
				SG->write_grouped_segments_shapefile(params->output_shapefile);
				SG->reconstruct_polygons(params->temp_dir);
			}

		}
	}
}