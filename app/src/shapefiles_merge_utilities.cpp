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

			template <typename point_type>
			bool pts_collinear(point_type p1, point_type p2, point_type p3) {
				double numerator = p2.y() * (p1.x() - p3.x()) + p1.y() * (p3.x() - p2.x()) + p3.y() * (p2.x() - p1.x());
				double denominator = (p2.x() - p1.x()) * (p1.x() - p3.x()) + (p2.y() - p1.y()) * (p1.y() - p3.y());
				double ratio = numerator / denominator;
				double points_angle = std::abs(std::fmod(std::atan(ratio), M_PI));
				points_angle = std::min<double>(points_angle, std::abs<double>(points_angle - M_PI));
				bool pts_collinear = (points_angle < (M_PI / 180));
				return pts_collinear;
			}

			double segments_overlap_ratios(const Inexact_Segment_2& segment1, const Inexact_Segment_2& segment2)
			{

				double seg1_Xs_min = std::min(segment1.source().x(), segment1.target().x());
				double seg2_Xs_min = std::min(segment2.source().x(), segment2.target().x());
				double seg1_Xs_max = std::max(segment1.source().x(), segment1.target().x());
				double seg2_Xs_max = std::max(segment2.source().x(), segment2.target().x());

				double seg1_Ys_min = std::min(segment1.source().y(), segment1.target().y());
				double seg2_Ys_min = std::min(segment2.source().y(), segment2.target().y());
				double seg1_Ys_max = std::max(segment1.source().y(), segment1.target().y());
				double seg2_Ys_max = std::max(segment2.source().y(), segment2.target().y());

				double Ax = std::min(seg1_Xs_max, seg2_Xs_max);
				double Ay = std::min(seg1_Ys_max, seg2_Ys_max);

				double Bx = std::max(seg1_Xs_min, seg2_Xs_min);
				double By = std::max(seg1_Ys_min, seg2_Ys_min);

				double X_overlap = std::max(double(0), Ax - Bx);
				double Y_overlap = std::max(double(0), Ay - By);

				auto safe_div = [](auto A, auto B) { return (B == 0) ? (1) : (A / B); };

				double line1_X_overlap = safe_div(X_overlap, (seg2_Xs_max - seg2_Xs_min));
				double line1_Y_overlap = safe_div(Y_overlap, (seg2_Ys_max - seg2_Ys_min));
				double line2_X_overlap = safe_div(X_overlap, (seg1_Xs_max - seg1_Xs_min));
				double line2_Y_overlap = safe_div(Y_overlap, (seg1_Ys_max - seg1_Ys_min));

				return  std::max(line2_X_overlap + line2_Y_overlap, line1_X_overlap + line1_Y_overlap);
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
									R.shrink_to_fit();

									//simplify edge
									std::vector<Inexact_Point_2> simplified_R;
									simplify_ring(R, simplified_R);

									for (size_t l = 0; l < simplified_R.size()-1; ++l) {
										all_segments.push_back(Inexact_Segment_2(simplified_R[l], simplified_R[l+1]));
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

			OGRLayer* create_output_layer(GDALDriver* driver,
				OGRSpatialReference* target_ref, boost::filesystem::path out_path,
				std::map<std::string, GDALDataset*>& temp_dataset_map) {

				GDALDataset* out_dataset = driver->Create(out_path.string().c_str(), 0, 0, 0, GDT_Unknown, NULL);
				temp_dataset_map[out_path.string()] = out_dataset;

				if (out_dataset == NULL) {
					throw std::ios_base::failure(fmt::format("Error : unable to write temporary shapefile at {}", out_path.string()));
				}

				OGRLayer* target_layer = out_dataset->CreateLayer(std::string("").c_str(), target_ref, wkbPolygon, NULL);
				if (target_layer == NULL) {
					throw std::logic_error("Error : layer creation failed.");
				}

				return target_layer;
			}

			void load_datasets_layers(std::vector<boost::filesystem::path> input_paths,
				std::map<size_t, GDALDataset*>& datasets_map,
				std::map<size_t, OGRLayer*>& layers_map) {
				

				for (size_t layer_index = 0; layer_index < input_paths.size(); ++layer_index)
				{
					datasets_map[layer_index] = (GDALDataset*)GDALOpenEx(input_paths[layer_index].string().c_str(), GDAL_OF_VECTOR, NULL, NULL, NULL);
				

				if (datasets_map[layer_index] == NULL) {
					throw std::ios_base::failure("Error : unable to load shapefile.");
				}

				layers_map[layer_index] = datasets_map[layer_index]->GetLayer(0);
				}
				
			}

			void clean_invalid(OGRLayer* c_layer) {

				size_t features_count = size_t(c_layer->GetFeatureCount());

				for (size_t c_feat_index = 0; c_feat_index < features_count; ++c_feat_index) {
					OGRFeature* feat = c_layer->GetNextFeature();
					if (feat == NULL) {
						c_layer->DeleteFeature(c_feat_index);
						continue;
					}

					OGRGeometry* geom = feat->GetGeometryRef();
					OGRwkbGeometryType geom_type = geom->getGeometryType();
					if (geom_type == wkbPolygon){
						continue;
					}
					else if (geom_type == wkbMultiPolygon) {
						OGRMultiPolygon* c_multipolygon = dynamic_cast<OGRMultiPolygon*>(geom);
						for (size_t sub_geom_idx = 0; sub_geom_idx < c_multipolygon->getNumGeometries(); ++sub_geom_idx) {
							OGRGeometry* sub_geom = c_multipolygon->getGeometryRef(sub_geom_idx);
							if (sub_geom->getGeometryType() == wkbPolygon) {
								OGRFeature* feature = OGRFeature::CreateFeature(c_layer->GetLayerDefn());
								feature->SetGeometry(sub_geom->clone());
								OGRErr error = c_layer->CreateFeature(feature);
								if (error != OGRERR_NONE) BOOST_LOG_TRIVIAL(warning) << fmt::format("Error code : {}", error);
								OGRFeature::DestroyFeature(feature);
							}
						}
						c_layer->DeleteFeature(c_feat_index);
						continue;
					}
					// delete all different from polygon and multipolygon
					else c_layer->DeleteFeature(c_feat_index);
				
				}
			}

			void simplify_ring(std::vector<Inexact_Point_2>& R, std::vector<Inexact_Point_2>& simplified_R) {
				simplified_R.reserve(R.size());
				Inexact_Point_2 turn_m2 = R[R.size() - 2];
				Inexact_Point_2 turn_m1 = R[R.size() - 1];
				for (size_t c_point_idx = 0; c_point_idx < R.size(); ++c_point_idx) {
					if (!pts_collinear(turn_m1, R[c_point_idx], turn_m2)) {
						simplified_R.push_back(turn_m1);
					}
					turn_m2 = turn_m1;
					turn_m1 = R[c_point_idx];
				}
				simplified_R.push_back(simplified_R[0]);
				simplified_R.shrink_to_fit();
			}

			void overlay_union_layers(std::vector<boost::filesystem::path>& regularized_layers_path) {

				//prepare parameters
				char** union_params = NULL;
				union_params = CSLAddString(union_params, "INPUT_PREFIX=o_");
				union_params = CSLAddString(union_params, "METHOD_PREFIX=n_");
				union_params = CSLAddString(union_params, "SKIP_FAILURES=YES");
				union_params = CSLAddString(union_params, "PROMOTE_TO_MULTI=NO");
				union_params = CSLAddString(union_params, "KEEP_LOWER_DIMENSION_GEOMETRIES=NO");
				std::map<std::string, GDALDataset*> temp_dataset_map;
				std::map<size_t, GDALDataset*> datasets_map;
				std::map<size_t, OGRLayer*> layers_map;

				try {
					load_datasets_layers(regularized_layers_path, datasets_map, layers_map);
				}
				catch (std::exception& e) {
					for (auto out_dataset = datasets_map.begin(); out_dataset != datasets_map.end(); ++out_dataset) {
						if (out_dataset->second != NULL) GDALClose(out_dataset->second);
					}
					BOOST_LOG_TRIVIAL(debug) << e.what();
					BOOST_LOG_TRIVIAL(fatal) << "Fatal error! For a quick resolution keep a copy of input data!";
					return;
				}

				try {
					const std::string driver_name = "ESRI Shapefile";

					GDALDriver* driver = GetGDALDriverManager()->GetDriverByName(driver_name.c_str());
					if (driver == NULL) {
						throw std::logic_error("Error : ESRI Shapefile driver not available.");
					}

					OGRSpatialReference* target_ref = layers_map[0]->GetSpatialRef();

					//// file_path constructions
					//boost::filesystem::path c_out_basename(fmt::format("union_{}.shp", regularized_layers_path.size()));
					//boost::filesystem::path outdir(params->temp_dir);
					//boost::filesystem::path c_full_path = outdir / c_out_basename;
					//OGRLayer* first_union_layer = create_output_layer(driver, target_ref, c_full_path, temp_dataset_map);
					//// union of first two layers
					//OGRErr err = layers_map[0]->Union(layers_map[1], first_union_layer, union_params, NULL, NULL);

					//if (err != OGRERR_NONE)
					//{
					//	BOOST_LOG_TRIVIAL(fatal) << "Union error of regularized layers";
					//	CSLDestroy(union_params);
					//	throw std::logic_error("OGR error at union step!");
					//}
					//first_union_layer->SyncToDisk();

					OGRLayer* last_union_layer = layers_map[0];
					for (size_t layer_index = 1; layer_index < regularized_layers_path.size(); ++layer_index) {
						boost::filesystem::path c_out_basename(fmt::format("union_{}.shp", regularized_layers_path.size()- layer_index));
						boost::filesystem::path outdir(params->temp_dir);
						boost::filesystem::path c_full_path = outdir / c_out_basename;
						OGRLayer* kth_union_layer = create_output_layer(driver, target_ref, c_full_path, temp_dataset_map);

						//union of kth layers
						OGRErr err = last_union_layer->Union(layers_map[layer_index], kth_union_layer, union_params, NULL, NULL);
						last_union_layer = kth_union_layer;

						if (err != OGRERR_NONE)
						{
							BOOST_LOG_TRIVIAL(fatal) << "Union error of regularized layers";
							CSLDestroy(union_params);
							throw std::logic_error("OGR error at union step!");
						}
						clean_invalid(last_union_layer);
						last_union_layer->SyncToDisk();
					}

				}
				catch (std::exception& e) {
					for (auto out_dataset = datasets_map.begin(); out_dataset != datasets_map.end(); ++out_dataset) {
						if (out_dataset->second != NULL) GDALClose(out_dataset->second);
					}
					for (auto out_dataset = temp_dataset_map.begin(); out_dataset != temp_dataset_map.end(); ++out_dataset) {
						if (out_dataset->second != NULL) GDALClose(out_dataset->second);
					}
					BOOST_LOG_TRIVIAL(debug) << e.what();
					BOOST_LOG_TRIVIAL(fatal) << "Fatal error! For a quick resolution keep a copy of input data!";
					return;
				}

				for (auto out_dataset = datasets_map.begin(); out_dataset != datasets_map.end(); ++out_dataset) {
					if (out_dataset->second != NULL) GDALClose(out_dataset->second);
				}
				for (auto out_dataset = temp_dataset_map.begin(); out_dataset != temp_dataset_map.end(); ++out_dataset) {
					if (out_dataset->second != NULL) GDALClose(out_dataset->second);
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