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
#include "geometry_lab.h"

namespace LxGeo
{
	namespace LxShapefilesMerge
	{
		namespace shapefiles_merge_utils
		{

			IK_to_EK to_exact;
			EK_to_IK to_inexact;

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

			

			short int cyclcic_order_distance(short int a, short int b, short int cycle) {
				short int result = abs(a - b)%cycle;
				return std::min<short int>(result, cycle - result);
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

			Point_2 exact_load(const Inexact_Point_2& in_pt) {
				std::stringstream stream;
				stream << std::setprecision(17) << in_pt.x() << " " << in_pt.y();
				FT x, y;
				stream >> x >> y;
				Point_2 out_pt = Point_2(x, y);
				return out_pt;
			}

			void load_segments_data(std::vector<std::string>& all_paths,
				std::vector<Segment_2>& all_segments,
				std::vector<short int>& segment_LID,
				std::vector<short int>& segment_PID,
				std::vector<short int>& segment_ORDinP,
				std::vector<short int>& segment_RRSize,
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

						tqdm bar;
						for (size_t j = 0; j < n; ++j) {
							bar.progress(j, n);
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
									std::vector<Inexact_Point_2> simplified_R = simplify_aberrant_ring(R.begin(), R.end());
									
									if (simplified_R.size() < 3) {
										//BOOST_LOG_TRIVIAL(debug) << "oversimplify " << i << "" << j;
										continue;
									}
									/*std::vector<Inexact_Point_2> simplified_R;
									simplify_ring(R, simplified_R);*/

									size_t c_ignored_segments_count = 0;
									for (size_t l = 0; l < simplified_R.size()-1; ++l) {
										
										Segment_2 segment_to_add(exact_load(simplified_R[l]), exact_load(simplified_R[l + 1]));
										if (segment_to_add.squared_length() == FT(0)) {
											BOOST_LOG_TRIVIAL(info) << "Ignoring segment having 0 length!";
											c_ignored_segments_count++;
											continue;
										}
										all_segments.push_back(segment_to_add);
										segment_LID.push_back(i);
										segment_PID.push_back(j);
										// if two segments share the same ORDinP & PID than at least an outer ring exists.
										segment_ORDinP.push_back(l- c_ignored_segments_count);
									}
									for (size_t added_segment_cnt = 0; added_segment_cnt < simplified_R.size() -1- c_ignored_segments_count; ++added_segment_cnt) {
										segment_RRSize.push_back(simplified_R.size() - 1 - c_ignored_segments_count);
									}

								}

							}
						}
						bar.finish();
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


			void make_rtree(const std::vector< Segment_2 >& all_segments, Boost_RTree_2& RT)
			{
				for (size_t i = 0; i < all_segments.size(); ++i) {
					const Segment_2 c_segment = all_segments[i];
										
					double xmin = std::min<double>(CGAL::to_double(c_segment.source().x()), CGAL::to_double(c_segment.target().x()));
					double xmax = std::max<double>(CGAL::to_double(c_segment.source().x()), CGAL::to_double(c_segment.target().x()));
					double ymin = std::min<double>(CGAL::to_double(c_segment.source().y()), CGAL::to_double(c_segment.target().y()));
					double ymax = std::max<double>(CGAL::to_double(c_segment.source().y()), CGAL::to_double(c_segment.target().y()));

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

			std::list<OGRPolygon> apply_fix_polygon(OGRPolygon* to_fix_polygon) {

				std::vector<Inexact_Point_2> exterior_ring_points = container_transform_OGRRING2vector_Points(to_fix_polygon->getExteriorRing());
				// could be a multitude of rings if explosion happened
				std::list<std::vector<Inexact_Point_2>> exterior_ring_fixed_pts = explose_self_intersecting_polygon(simplify_aberrant_ring(exterior_ring_points.begin(), exterior_ring_points.end()));
				std::list<std::vector<Inexact_Point_2>> interior_rings_fixed_pts;
				if (exterior_ring_fixed_pts.size() == 1) {
					for (int k = 0; k < to_fix_polygon->getNumInteriorRings(); ++k) {
						auto int_ring = container_transform_OGRRING2vector_Points(to_fix_polygon->getInteriorRing(k));
						interior_rings_fixed_pts.push_back(
							simplify_aberrant_ring(int_ring.begin(), int_ring.end())
						);
					}
				}

				std::list<OGRPolygon> fixed_polygons_list;
				for (std::vector<Inexact_Point_2> exploded_ring_pts : exterior_ring_fixed_pts) {					
					OGRPolygon current_polygon;
					OGRLinearRing copy_ring = container_transform_vector_Points2OGRRING(exploded_ring_pts);
					current_polygon.addRing(&copy_ring);
					for (std::vector<Inexact_Point_2> c_interior_ring : interior_rings_fixed_pts) {
						OGRLinearRing copy_interior_ring = container_transform_vector_Points2OGRRING(c_interior_ring);
						current_polygon.addRing(&copy_interior_ring);
					}
					fixed_polygons_list.push_back(current_polygon);			
				}				

				return fixed_polygons_list;
			}

			std::list<OGRPolygon> apply_fix_polygon_temp(OGRPolygon* to_fix_polygon) {

				OGRPolygon current_polygon;

				std::vector<Inexact_Point_2> exterior_ring_points = container_transform_OGRRING2vector_Points(to_fix_polygon->getExteriorRing());
				std::vector<Inexact_Point_2> exterior_ring_fixed_pts = simplify_aberrant_ring(exterior_ring_points.begin(), exterior_ring_points.end());
				OGRLinearRing copy_ring = container_transform_vector_Points2OGRRING(exterior_ring_points);
				current_polygon.addRing(&copy_ring);

				std::list<std::vector<Inexact_Point_2>> interior_rings_fixed_pts;
				for (int k = 0; k < to_fix_polygon->getNumInteriorRings(); ++k) {
					auto inner_ring = container_transform_OGRRING2vector_Points(to_fix_polygon->getInteriorRing(k));
					interior_rings_fixed_pts.push_back(
						simplify_aberrant_ring(inner_ring.begin(), inner_ring.end())
					);
				}

				for (std::vector<Inexact_Point_2> c_interior_ring : interior_rings_fixed_pts) {
					OGRLinearRing copy_interior_ring = container_transform_vector_Points2OGRRING(c_interior_ring);
					current_polygon.addRing(&copy_interior_ring);
				}

				std::list<OGRPolygon> fixed_polygons_list;
				fixed_polygons_list.push_back(current_polygon);
				return fixed_polygons_list;
			}

			void fix_layer_invalid_geometries(OGRLayer* c_layer) {

				size_t features_count = size_t(c_layer->GetFeatureCount());

				std::list<OGRPolygon> polygons_to_add;
				std::list<size_t> features_to_remove;

				for (size_t c_feat_index = 0; c_feat_index < features_count; ++c_feat_index) {
					OGRFeature* feat = c_layer->GetNextFeature();
					if (feat == NULL) {
						//c_layer->DeleteFeature(c_feat_index);
						features_to_remove.push_back(c_feat_index);
						continue;
					}

					OGRGeometry* geom = feat->GetGeometryRef();
					OGRwkbGeometryType geom_type = geom->getGeometryType();
					if (geom_type == wkbPolygon){
						OGRPolygon* P = dynamic_cast<OGRPolygon*>(geom);
						std::list<OGRPolygon> fixed_Ps = apply_fix_polygon_temp(P);
						polygons_to_add.splice(polygons_to_add.end(), fixed_Ps);
						features_to_remove.push_back(c_feat_index);
						continue;
					}
					else if (geom_type == wkbMultiPolygon) {
						OGRMultiPolygon* c_multipolygon = dynamic_cast<OGRMultiPolygon*>(geom);
						for (size_t sub_geom_idx = 0; sub_geom_idx < c_multipolygon->getNumGeometries(); ++sub_geom_idx) {
							OGRGeometry* sub_geom = c_multipolygon->getGeometryRef(sub_geom_idx);
							if (sub_geom->getGeometryType() == wkbPolygon) {
								// get current polygon rings
								OGRPolygon* P = dynamic_cast<OGRPolygon*>(sub_geom);
								std::list<OGRPolygon> fixed_Ps = apply_fix_polygon_temp(P);
								polygons_to_add.splice(polygons_to_add.end(), fixed_Ps);

							}
						}
						features_to_remove.push_back(c_feat_index);
						continue;
					}
					// delete all different from polygon and multipolygon
					else features_to_remove.push_back(c_feat_index);
				
				}

				// remove all selected features
				for (size_t f_to_remove : features_to_remove) {
					c_layer->DeleteFeature(f_to_remove);
				}
				// add all fixed geometries
				for (OGRPolygon poly_to_add : polygons_to_add) {
					OGRFeature* feature = OGRFeature::CreateFeature(c_layer->GetLayerDefn());
					feature->SetGeometry(poly_to_add.clone());
					OGRErr error = c_layer->CreateFeature(feature);
					if (error != OGRERR_NONE) BOOST_LOG_TRIVIAL(warning) << fmt::format("Error code : {}", error);
					OGRFeature::DestroyFeature(feature);
				}

				c_layer->ResetReading();
			}



			std::string overlay_union_layers(std::vector<boost::filesystem::path>& regularized_layers_path) {

				std::string last_overlay_layer_path;
				//prepare parameters
				char** union_params = NULL;
				union_params = CSLAddString(union_params, "INPUT_PREFIX=o_");
				union_params = CSLAddString(union_params, "METHOD_PREFIX=n_");
				union_params = CSLAddString(union_params, "SKIP_FAILURES=YES");
				union_params = CSLAddString(union_params, "USE_PREPARED_GEOMETRIES=NO");
				union_params = CSLAddString(union_params, "KEEP_LOWER_DIMENSION_GEOMETRIES=YES");
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
					return "";
				}

				try {
					const std::string driver_name = "ESRI Shapefile";

					GDALDriver* driver = GetGDALDriverManager()->GetDriverByName(driver_name.c_str());
					if (driver == NULL) {
						throw std::logic_error("Error : ESRI Shapefile driver not available.");
					}

					OGRSpatialReference* target_ref = layers_map[0]->GetSpatialRef();

					OGRLayer* last_union_layer = layers_map[0];
					for (size_t layer_index = 1; layer_index < regularized_layers_path.size(); ++layer_index) {
						boost::filesystem::path c_out_basename(fmt::format("union_{}.shp", regularized_layers_path.size()- layer_index));
						boost::filesystem::path outdir(params->temp_dir);
						boost::filesystem::path c_full_path = outdir / c_out_basename;
						last_overlay_layer_path = c_full_path.string();
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
						fix_layer_invalid_geometries(last_union_layer);
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
					return "";
				}

				for (auto out_dataset = datasets_map.begin(); out_dataset != datasets_map.end(); ++out_dataset) {
					if (out_dataset->second != NULL) GDALClose(out_dataset->second);
				}
				for (auto out_dataset = temp_dataset_map.begin(); out_dataset != temp_dataset_map.end(); ++out_dataset) {
					if (out_dataset->second != NULL) GDALClose(out_dataset->second);
				}
				
				return last_overlay_layer_path;
			}


			void regularize_segments(std::vector<Segment_2>& all_segments,
				std::vector<short int>& segment_LID,
				std::vector<short int>& segment_PID,
				std::vector<short int>& segment_ORDinP,
				std::vector<short int>& segment_RRSize) {

				BOOST_LOG_TRIVIAL(info) << "Computing segments angles!";
				std::vector<double> segment_angle;
				segment_angle.reserve(all_segments.size());
				for (auto& c_segment : all_segments) {
					Inexact_Point_2 p1(CGAL::to_double(c_segment.source().x()) , CGAL::to_double((c_segment.source().y())));
					Inexact_Point_2 p2(CGAL::to_double(c_segment.target().x()), CGAL::to_double(c_segment.target().y()));
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
					segment_LID, segment_PID, segment_ORDinP, segment_RRSize, segment_angle, segments_tree);

				//SG->write_grouped_segments_shapefile(params->output_shapefile);
				SG->fill_graph();
				SG->cluster_segments();
				SG->fuse_segments();
				SG->write_grouped_segments_shapefile(params->output_shapefile);
				SG->reconstruct_polygons(params->temp_dir);

			}


			void read_shapefiles(const std::vector<std::string>& all_paths, std::vector<std::vector<std::vector<std::vector<Inexact_Point_2> > > >& all_polygons)
			{
				if (all_paths.empty()) return;

				GDALDataset* dataset_ref = NULL;
				GDALDataset* dataset = NULL;

				try {
					// Gets a spatial reference
					dataset_ref = (GDALDataset*)GDALOpenEx(all_paths[0].c_str(), GDAL_OF_VECTOR | GDAL_OF_UPDATE, NULL, NULL, NULL);
					if (dataset_ref == NULL) {
						throw std::ios_base::failure("Error : unable to load shapefile.");
					}
					OGRLayer* layer_ref = dataset_ref->GetLayer(0);
					OGRSpatialReference* target_ref = layer_ref->GetSpatialRef();

					all_polygons.reserve(all_paths.size());
					for (size_t i = 0; i < all_paths.size(); ++i) {
						std::vector<std::vector<std::vector<Inexact_Point_2> > > polygons;
						if (i == 0) {
							read_single_shapefile(dataset_ref, target_ref, polygons);
						}
						else {
							read_single_shapefile(all_paths[i], target_ref, polygons);
						}
						
						all_polygons.push_back(polygons);
					}
				}
				catch (std::exception& e) {
					std::cerr << e.what() << std::endl;
					all_polygons.clear();
				}

				if (dataset_ref != NULL) GDALClose(dataset_ref);
				if (dataset != NULL) GDALClose(dataset);
			}

			void read_single_shapefile(std::string shapefile_path, OGRSpatialReference* target_ref, std::vector<std::vector<std::vector<Inexact_Point_2> > >& polygons) {
				GDALDataset* dataset = NULL;
				dataset = (GDALDataset*)GDALOpenEx(shapefile_path.c_str(), GDAL_OF_VECTOR | GDAL_OF_UPDATE, NULL, NULL, NULL);
				read_single_shapefile(dataset, target_ref, polygons);
				// Reading is over, closes the dataset
				GDALClose(dataset);
			}

			void read_single_shapefile(GDALDataset* shapefile_dataset, OGRSpatialReference* target_ref, std::vector<std::vector<std::vector<Inexact_Point_2> > >& polygons) {
				GDALDataset* dataset = NULL;

				dataset = shapefile_dataset;				

				if (dataset == NULL) {
					throw std::ios_base::failure("Error : unable to load shapefile.");
				}

				OGRLayer* layer = dataset->GetLayer(0);
				fix_layer_invalid_geometries(layer);
				layer->SyncToDisk();
				OGRSpatialReference* source_ref = layer->GetSpatialRef();
				OGRCoordinateTransformation* CT = OGRCreateCoordinateTransformation(source_ref, (target_ref==NULL)? source_ref : target_ref);

				// Step 2.
				// Reads contents of the current shapefile

				size_t n = size_t(layer->GetFeatureCount());
				polygons.reserve(n);

				for (size_t j = 0; j < n; ++j) {
					OGRFeature* feat = layer->GetNextFeature();
					if (feat == NULL) continue;

					OGRGeometry* geom = feat->GetGeometryRef();

					// Assumes the shapefile only contains OGRPolygons
					if (OGRPolygon* P = dynamic_cast<OGRPolygon*>(geom)) {
						std::vector<std::vector<Inexact_Point_2> > polygon;

						// Reads rings
						std::list<OGRLinearRing*> ogr_rings;
						ogr_rings.push_back(P->getExteriorRing());
						for (int k = 0; k < P->getNumInteriorRings(); ++k) ogr_rings.push_back(P->getInteriorRing(k));

						for (OGRLinearRing* ring : ogr_rings) {
							std::vector<Inexact_Point_2> R;
							int ring_size = ring->getNumPoints();

							R.reserve(ring_size);
							for (int u = 0; u < ring_size ; ++u) { //-1
								OGRPoint pt;
								ring->getPoint(u, &pt);
								double x = pt.getX(), y = pt.getY();
								CT->Transform(1, &x, &y);
								R.push_back(Inexact_Point_2(x, y));
							}

							polygon.push_back(R);
						}

						polygons.push_back(polygon);
					}
				}

			}

			void extract_edges_from_polygons(std::vector<std::vector<std::vector<Inexact_Point_2>>>& c_layer_regularized_polygons,
				std::vector<Boost_LineString_2>& c_layer_regularized_edges) {

				std::list<Boost_LineString_2> linestrings_list;

				for (auto c_polygon : c_layer_regularized_polygons) {
					for (auto c_ring : c_polygon) {
						
						for (size_t c_pt_idx = 0; c_pt_idx < c_ring.size(); ++c_pt_idx) {
							size_t next_pt_idx = (c_pt_idx + 1) % c_ring.size();
							Inexact_Point_2 src_pt = c_ring[c_pt_idx];
							Inexact_Point_2 trg_pt = c_ring[next_pt_idx];
							Boost_LineString_2 c_linestring;
							bg::append(c_linestring, Boost_Point_2(src_pt.x(), src_pt.y()));
							bg::append(c_linestring, Boost_Point_2(trg_pt.x(), trg_pt.y()));
							linestrings_list.push_back(c_linestring);
						}
						
					}
				}
				c_layer_regularized_edges = std::vector<Boost_LineString_2>(linestrings_list.cbegin(), linestrings_list.cend());
			}

			void make_rtree_polygons(const std::vector<std::vector<std::vector<Inexact_Point_2> > >& polygons, Boost_RTree_2& RT)
			{
				for (size_t i = 0; i < polygons.size(); ++i) {
					const std::vector<std::vector<Inexact_Point_2> >& polygon = polygons[i];

					// Computes bbox of outer ring

					Bbox_2 B = CGAL::bbox_2(polygon[0].cbegin(), polygon[0].cend());
					double xmin = B.xmin(), xmax = B.xmax(), ymin = B.ymin(), ymax = B.ymax();

					Boost_Box_2 box(Boost_Point_2(xmin, ymin), Boost_Point_2(xmax, ymax));
					RT.insert(Boost_Value_2(box, i));
				}
			}

			void make_rtree_linestrings(const std::vector<Boost_LineString_2>& linestrings, Boost_RTree_2& RT)
			{
				for (size_t i = 0; i < linestrings.size(); ++i) {
					const Boost_LineString_2& c_linestring = linestrings[i];
					Boost_Box_2 envelope;
					bg::envelope(c_linestring, envelope);
					RT.insert(Boost_Value_2(envelope, i));
				}
			}

		}
	}
}