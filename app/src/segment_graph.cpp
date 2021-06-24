#include "segment_graph.h"
#include "defs.h"
#include "defs_cgal.h"
#include "defs_boost.h"
#include "parameters.h"
#include "segment_graph_wrapper.h"
#include <boost/graph/graphviz.hpp>
#include <boost/filesystem.hpp>
#include <fmt/core.h>
#include <CGAL/linear_least_squares_fitting_2.h>
#include <assert.h>
#include <limits>
#include <gdal.h>
#include <gdal_priv.h>
#include <ogrsf_frmts.h>
#include <CGAL/exceptions.h>


namespace LxGeo
{
	namespace LxShapefilesMerge
	{

		
		SegmentGraph::SegmentGraph(std::vector<Inexact_Segment_2>& all_segments,
			std::vector<short int>& segment_LID,
			std::vector<short int>& segment_PID,
			std::vector<short int>& segment_ORDinP,
			std::vector<double>& segment_angle,
			Boost_RTree_2& segments_tree): _all_segments(all_segments), _segment_LID(segment_LID), _segment_PID(segment_PID), _segment_ORDinP(segment_ORDinP), _segment_angle(segment_angle),
			_segments_tree(segments_tree)
		{
			// load weights and thresh from parameters
			e_distance_weight = params->e_distance_weight;
			e_angle_weight = params->e_angle_weight;
			MAX_GROUPING_DISTANCE = params->MAX_GROUPING_DISTANCE;
			MAX_GROUPING_ANGLE_DEG = params->MAX_GROUPING_ANGLE_DEG;
			MAX_GROUPING_ANGLE_RAD = MAX_GROUPING_ANGLE_DEG * M_PI / 180;

			
			SG = BoostSegmentGraph(all_segments.size());
			vertcies_centrality = std::vector<double>(all_segments.size());
			vertcies_groups = std::vector<size_t>(all_segments.size(),0);
			groupes_count = 1;
		}


		SegmentGraph::~SegmentGraph()
		{
		}

		void SegmentGraph::fill_graph()
		{

			vertcies_centrality = std::vector<double>(_all_segments.size());


			for (size_t c_segment_index=0; c_segment_index <_all_segments.size(); ++c_segment_index)
			{

				std::vector<size_t> neighborhood_indices;
				std::vector<double> neighborhood_weights;
				std::vector<double> neighborhood_distance;

				get_neighborhood_segments_indices_weights(c_segment_index,
					neighborhood_indices,
					neighborhood_weights,
					neighborhood_distance
				);

				//test only
				assert (neighborhood_indices.size() == neighborhood_weights.size());

				double sum_vertex_centrality = 0;
				// adding graph edges
				for (size_t neigh_iter_index=0; neigh_iter_index < neighborhood_indices.size(); ++neigh_iter_index)
				{

					size_t neigh_index = neighborhood_indices[neigh_iter_index];
					double neigh_weight = neighborhood_weights[neigh_iter_index];
					double neigh_distance = neighborhood_distance[neigh_iter_index];
					

					auto add_pair = boost::add_edge(vertex_descriptor(c_segment_index), vertex_descriptor(neigh_index), SG);
					auto e = add_pair.first;
					//SG[e].distance = neigh_distance;

					// the weight can be assigned using a weight map.
					/*boost::property_map<BoostSegmentGraph, boost::edge_weight_t>::type weightmap =
						get(boost::edge_weight, SG);
					weightmap[e] = neigh_weight;*/
					sum_vertex_centrality += neigh_weight;

				}

				/*boost::property_map<BoostSegmentGraph, boost::vertex_index_t>::type vertex_index_map =
					boost::get(boost::vertex_index, SG);

				vertex_index_map[c_segment_index] = c_segment_index;

				boost::property_map<BoostSegmentGraph, boost::vertex_centrality_t>::type centrality_map =
					boost::get(boost::vertex_centrality, SG);
				*/

				double normalized_centrality_value = (sum_vertex_centrality / neighborhood_indices.size());
				//centrality_map[c_segment_index] = normalized_centrality_value;
				vertcies_centrality[c_segment_index]= normalized_centrality_value;
			}

		}

		void SegmentGraph::get_neighborhood_segments_indices_weights(size_t c_segment_index,
			std::vector<size_t>& neighborhood_indices,
			std::vector<double>& neighborhood_weights,
			std::vector<double>& neighborhood_distance
		)
		{
			Inexact_Segment_2 c_segment = _all_segments[c_segment_index];
			// First filter (distance search)
			std::vector<size_t> distance_neighborhood_indices;

			double x_min, x_max, y_min, y_max;
			double p1_x, p2_x, p1_y, p2_y;
			p1_x = c_segment.source().x();
			p2_x = c_segment.target().x();
			p1_y = c_segment.source().y();
			p2_y = c_segment.target().y();

			x_min = std::min<double>(p1_x, p2_x) - MAX_GROUPING_DISTANCE;
			x_max = std::max<double>(p1_x, p2_x) + MAX_GROUPING_DISTANCE;
			y_min = std::min<double>(p1_y, p2_y) - MAX_GROUPING_DISTANCE;
			y_max = std::max<double>(p1_y, p2_y) + MAX_GROUPING_DISTANCE;

			Boost_Segment_2 c_boost_segment = Boost_Segment_2(
				Boost_Point_2(p1_x, p1_y),
				Boost_Point_2(p2_x, p2_y)
			);

			std::vector<Boost_Value_2> possible_neighbors;
			Boost_Box_2 query(Boost_Point_2(x_min, y_min), Boost_Point_2(x_max, y_max));

			_segments_tree.query(bgi::intersects(query), std::back_inserter(possible_neighbors));

			for (int i = 0; i < possible_neighbors.size(); i++) {

				size_t possible_i_index = possible_neighbors[i].second;
				Inexact_Segment_2 possible_i_segment = _all_segments[possible_i_index];

				Boost_Segment_2 possible_i_boost_segment = Boost_Segment_2(
					Boost_Point_2(possible_i_segment.source().x(), possible_i_segment.source().y()),
					Boost_Point_2(possible_i_segment.target().x(), possible_i_segment.target().y())
				);
				//check distance
				double diff_distance = std::abs<double>(bg::distance(c_boost_segment, possible_i_boost_segment));
				if (diff_distance < MAX_GROUPING_DISTANCE)
				{
					// check angle difference
					double diff_angle = std::abs<double>(_segment_angle[possible_i_index] - _segment_angle[c_segment_index]);
					if (diff_angle < MAX_GROUPING_ANGLE_RAD)
					{
						//// ADD segment overlap ratios filter

						neighborhood_indices.push_back(possible_i_index);
						neighborhood_distance.push_back(diff_distance);
						double edge_pre_weight = ((diff_distance / MAX_GROUPING_DISTANCE) + (diff_angle / MAX_GROUPING_ANGLE_RAD)) / 2;
						// below is to decrease edge weight based on LID & PID of segments
						if (_segment_LID[c_segment_index] != _segment_LID[possible_i_index]) neighborhood_weights.push_back(edge_pre_weight / 3);
						else if (_segment_PID[c_segment_index] != _segment_PID[possible_i_index]) neighborhood_weights.push_back(edge_pre_weight / 2);
						else neighborhood_weights.push_back(edge_pre_weight);
					}
				}
			}
		}

		void SegmentGraph::cluster_segments()
		{
			
			size_t grouped_vertcies_count = 0;

			
			while (grouped_vertcies_count < _all_segments.size()) {
				// Search for minimum centrality_index as a starting group central segment
				
				size_t c_central_vertex_node = get_min_centrality_vertex_index();
				std::vector<size_t> c_connected_vertices_indices;
				get_connected_vertices_indices(c_central_vertex_node, c_connected_vertices_indices);

				size_t newly_added_count = create_segment_group(groupes_count, c_connected_vertices_indices);

				grouped_vertcies_count += newly_added_count;
				groupes_count++;

			}


		}

		size_t SegmentGraph::get_min_centrality_vertex_index()
		{
			return std::min_element(vertcies_centrality.begin(), vertcies_centrality.end()) - vertcies_centrality.begin();
		}

		void SegmentGraph::get_connected_vertices_indices(size_t vertex_idx, std::vector<size_t>& c_connected_vertices_indices) {

			boost::property_map<BoostSegmentGraph, boost::vertex_index_t>::type vertex_index_map =
				boost::get(boost::vertex_index, SG);

			c_connected_vertices_indices.push_back( vertex_idx );
			adjacency_iterator ai, ai_end;
			for (tie(ai, ai_end) = boost::adjacent_vertices(vertex_descriptor(vertex_idx), SG); ai != ai_end; ++ai) {
				c_connected_vertices_indices.push_back(get(vertex_index_map, *ai));
			}

				
		}

		size_t SegmentGraph::create_segment_group(size_t c_group_id, std::vector<size_t> c_connected_vertices_indices)
		{
			groups_map[c_group_id] = std::vector<size_t>();
			groups_map[c_group_id].reserve(c_connected_vertices_indices.size());
			size_t newly_added_count = 0;
			for (auto segment_id = c_connected_vertices_indices.begin(); segment_id != c_connected_vertices_indices.end(); ++segment_id)
			{
				if (vertcies_groups[*segment_id] == 0)
				{
					groups_map[c_group_id].push_back(*segment_id);
					vertcies_groups[*segment_id] = c_group_id;
					newly_added_count++;
				}

				boost::clear_vertex(vertex_descriptor(*segment_id), SG);
				delete_vertex_related(*segment_id);
			}
			groups_map[c_group_id].shrink_to_fit();

			return newly_added_count;
		}

		void SegmentGraph::delete_vertex_related(size_t vertex_idx)
		{
			//Assign high value for vertex_centrality
			vertcies_centrality[vertex_idx] = DBL_MAX;

			//remove vertex from graph
			//boost::remove_vertex(vertex_idx, SG);

		}

		void SegmentGraph::write_grouped_segments_shapefile(const std::string& output_filename)
		{
			if (_all_segments.empty()) {
				std::cout << "Warning : empty list of segments. No output written." << std::endl;
				return;
			}

			GDALDataset* source_dataset = NULL;
			GDALDataset* target_dataset = NULL;

			try {
				const std::string driver_name = "ESRI Shapefile";

				GDALDriver* driver = GetGDALDriverManager()->GetDriverByName(driver_name.c_str());
				if (driver == NULL) {
					throw std::logic_error("Error : ESRI Shapefile driver not available.");
				}

				// Step 1.
				// Reopens source file to get features

				source_dataset = (GDALDataset*)GDALOpenEx(params->paths[0].c_str(), GDAL_OF_VECTOR, NULL, NULL, NULL);
				if (source_dataset == NULL) {
					throw std::logic_error("Error : cannot reopen source shapefile.");
				}

				OGRLayer* source_layer = source_dataset->GetLayer(0);

				// Step 2.
				// Writes target file

				target_dataset = driver->Create(output_filename.c_str(), 0, 0, 0, GDT_Unknown, NULL);
				if (target_dataset == NULL) {
					throw std::logic_error("Error : creation of output file failed.");
				}

				OGRLayer* target_layer = target_dataset->CreateLayer(source_layer->GetName(), source_layer->GetSpatialRef(), wkbLineString, NULL);
				if (target_layer == NULL) {
					throw std::logic_error("Error : layer creation failed.");
				}

				OGRFieldDefn o_field_id("ID", OFTInteger);

				if (target_layer->CreateField(&o_field_id) != OGRERR_NONE) {
					throw std::logic_error("Error : field creation failed.");
				}

				OGRFieldDefn o_field_c_group_id("c_group_id", OFTInteger);

				if (target_layer->CreateField(&o_field_c_group_id) != OGRERR_NONE) {
					throw std::logic_error("Error : field c_group_id creation failed.");
				}

				OGRFieldDefn o_field_angle("angle", OFTReal);

				if (target_layer->CreateField(&o_field_angle) != OGRERR_NONE) {
					throw std::logic_error("Error : field angle creation failed.");
				}

				for (size_t i = 0; i < _all_segments.size(); ++i) {
					Inexact_Segment_2 S = _all_segments[i];
					OGRLineString ogr_linestring;
										
					const Inexact_Point_2& S1 = S.source();
					ogr_linestring.addPoint(&OGRPoint(S1.x(), S1.y()));
					const Inexact_Point_2& S2 = S.target();
					ogr_linestring.addPoint(&OGRPoint(S2.x(), S2.y()));
					
					OGRFeature* feature;
					feature = OGRFeature::CreateFeature(target_layer->GetLayerDefn());

					feature->SetGeometry(&ogr_linestring);
					feature->SetField("ID", int(i));
					feature->SetField("c_group_id", int(vertcies_groups[i]));
					feature->SetField("angle", _segment_angle[i]);

					// Writes new feature
					OGRErr error = target_layer->CreateFeature(feature);
					if (error != OGRERR_NONE) std::cout << "Error code : " << int(error) << std::endl;
					OGRFeature::DestroyFeature(feature);
				}

			}
			catch (std::exception& e) {
				std::cout << e.what() << std::endl;
			}

			if (source_dataset != NULL) GDALClose(source_dataset);
			if (target_dataset != NULL) GDALClose(target_dataset);
		}

		void SegmentGraph::fuse_segments()
		{
			for (size_t c_groupe_idx = 1; c_groupe_idx <= groupes_count; ++c_groupe_idx)
			{
				try {
					std::vector<Inexact_Segment_2*> respective_segments;
					get_respective_segment(c_groupe_idx, respective_segments);
				
					Inexact_Line_2 fitted_line;
					//get_best_fitting_line(fitted_line, respective_segments);
					double longest_value = 0;
					for (auto segment_i : respective_segments) {
						if (segment_i->squared_length() > longest_value) {
							fitted_line = segment_i->supporting_line();
							longest_value = segment_i->squared_length();
						}
					}

					for (size_t resp_idx = 0; resp_idx < respective_segments.size(); ++resp_idx) {

						Inexact_Segment_2 c_segment = *(respective_segments[resp_idx]);
						size_t idx_in_all_segs = groups_map[c_groupe_idx][resp_idx];
						Inexact_Point_2 source_projected = fitted_line.projection(c_segment.source());
						Inexact_Point_2 target_projected = fitted_line.projection(c_segment.target());
						// not sure about line below // maybe new Inexact_Segment_2 //fix func args

						_all_segments[idx_in_all_segs] = Inexact_Segment_2(source_projected, target_projected);
					}
				}
				catch (std::exception& e) {
					BOOST_LOG_TRIVIAL(debug) << "cgal error";
				}


			}
		}

		void SegmentGraph::get_respective_segment(size_t groupe_idx, std::vector<Inexact_Segment_2*>& respective_segments)
		{
			std::vector<size_t> respective_segments_indices = groups_map[groupe_idx];
			for (size_t idx = 0; idx < respective_segments_indices.size(); ++idx)
			{
				respective_segments.push_back(&(_all_segments[respective_segments_indices[idx]]));
			}
		}

		void SegmentGraph::get_best_fitting_line(Inexact_Line_2& fitted_line, std::vector<Inexact_Segment_2*>& respective_segments)
		{
			// init 1 weighted vector
			std::vector<size_t> segments_weights(respective_segments.size(), 1);
			if (respective_segments.size() == 1)
			{
				fitted_line = Inexact_Line_2(respective_segments[0]->source(), respective_segments[0]->source());
			}
			else
			{ 
				get_best_fitting_line(fitted_line, respective_segments, segments_weights);
			}

		}
		void SegmentGraph::get_best_fitting_line(Inexact_Line_2& fitted_line, std::vector<Inexact_Segment_2*>& respective_segments, std::vector<size_t>& segments_weights)
		{
			std::list<Inexact_Segment_2> all_segments;
			//test only
			assert(respective_segments.size() == segments_weights.size());
			// fill all_points with repetition (based on weights)
			for (size_t c_idx = 0; c_idx < respective_segments.size(); ++c_idx)
			{
				size_t c_segment_weight = segments_weights[c_idx];
				for (size_t rep = 0; rep < c_segment_weight; ++rep) {					
					all_segments.push_back(*respective_segments[c_idx]);
				}
			}
			// use cgal fit_line
			CGAL::linear_least_squares_fitting_2(all_segments.begin(), all_segments.end(), fitted_line, CGAL::Dimension_tag<1>());

		}

		void SegmentGraph::reconstruct_polygons(const std::string& temp_dir)
		{


			GDALDataset* source_dataset = NULL;
			std::map<size_t, GDALDataset*> outlayers_datasets_map;

			try {

				const std::string driver_name = "ESRI Shapefile";

				GDALDriver* driver = GetGDALDriverManager()->GetDriverByName(driver_name.c_str());
				if (driver == NULL) {
					throw std::logic_error("Error : ESRI Shapefile driver not available.");
				}

				// Gets a first spatial reference
				source_dataset = (GDALDataset*)GDALOpenEx(params->paths[0].c_str(), GDAL_OF_VECTOR, NULL, NULL, NULL);
				if (source_dataset == NULL) {
					throw std::ios_base::failure("Error : unable to load shapefile.");
				}
				OGRLayer* source_layer = source_dataset->GetLayer(0);
				OGRSpatialReference* target_ref = source_layer->GetSpatialRef();

				//Create output layers datasets			

				for (size_t layer_id = 0; layer_id < params->paths.size(); ++layer_id) {

					// file_path constructions
					boost::filesystem::path c_out_basename(fmt::format("l_{}.shp", layer_id));
					boost::filesystem::path outdir(params->temp_dir);
					boost::filesystem::path c_full_path = outdir / c_out_basename;

					outlayers_datasets_map[layer_id] = driver->Create(c_full_path.string().c_str(), 0, 0, 0, GDT_Unknown, NULL);
					if (outlayers_datasets_map[layer_id] == NULL) {
						throw std::ios_base::failure(fmt::format("Error : unable to write temporary shapefile at {}", c_full_path.string()));
					}

					OGRLayer* target_layer = outlayers_datasets_map[layer_id]->CreateLayer(source_layer->GetName(), target_ref, wkbPolygon, NULL);
					if (target_layer == NULL) {
						throw std::logic_error("Error : layer creation failed.");
					}

					OGRFieldDefn o_field_id("ID", OFTInteger);

					if (target_layer->CreateField(&o_field_id) != OGRERR_NONE) {
						throw std::logic_error("Error : field creation failed.");
					}
				}

					OGRLayer* current_dataset_layer= outlayers_datasets_map[0]->GetLayer(0);
					OGRFeature* feature;
					std::list<OGRLinearRing> ex_int_rings;
					OGRLinearRing current_ring;
					short int last_s_ORDinP=-1;
					short int last_s_PID=0;
					short int last_s_LID=0;

					for (size_t current_segment_index = 0; current_segment_index < (_all_segments.size() - 1); ++current_segment_index) {

						/*c_s_ORDinP = _segment_ORDinP[current_segment_index];
						c_s_PID = _segment_PID[current_segment_index];
						c_s_LID = _segment_LID[current_segment_index];*/

						if(_segment_ORDinP[current_segment_index] ==0)//changement of ring
						{
							if (!current_ring.IsEmpty()) {
								//closing ring by adding first element again
								/*OGRPoint first_point;
								current_ring.getPoint(0, &first_point);
								current_ring.addPoint(&OGRPoint(first_point));*/
								current_ring.closeRings();
								//add it to list of current polygon rings
								ex_int_rings.push_back(OGRLinearRing(current_ring));
								current_ring.empty();
							}


							if (_segment_LID[current_segment_index] != last_s_LID) {
								current_dataset_layer = outlayers_datasets_map[_segment_LID[current_segment_index]]->GetLayer(0);
								last_s_LID = _segment_LID[current_segment_index];
							}

							if (_segment_PID[current_segment_index] != last_s_PID) {

								OGRPolygon current_polygon;

								for (OGRLinearRing ring : ex_int_rings)
								{
									OGRLinearRing copy_ring(ring);
									current_polygon.addRing(&copy_ring);
								}
								//write polygon to respective dataset
								feature = OGRFeature::CreateFeature(current_dataset_layer->GetLayerDefn());
								feature->SetGeometry(&current_polygon);
								feature->SetField("ID", _segment_PID[current_segment_index]);
								// saving feature
								OGRErr error = current_dataset_layer->CreateFeature(feature);
								if (error != OGRERR_NONE) BOOST_LOG_TRIVIAL(warning)<< fmt::format("Error code : {}", error);
								OGRFeature::DestroyFeature(feature);
								// clearing old rings
								ex_int_rings.clear();
								last_s_PID = _segment_PID[current_segment_index];
							}

						}

						// adding point to rings (in all cases)
						Inexact_Point_2 s1_target = _all_segments[current_segment_index].target();
						Inexact_Segment_2 next_segment_in_polygon;
						// below is to assure next segment is in the same polygon
						if (_segment_ORDinP[current_segment_index + 1] != 0) next_segment_in_polygon = _all_segments[current_segment_index + 1];
						else next_segment_in_polygon = _all_segments[current_segment_index - _segment_ORDinP[current_segment_index]];

						Inexact_Point_2 projected_point = next_segment_in_polygon.supporting_line().projection(s1_target);
						current_ring.addPoint(&OGRPoint(projected_point.x(), projected_point.y()));
					}

				

			}
			catch (std::exception& e) {
				BOOST_LOG_TRIVIAL(debug) << e.what();
				BOOST_LOG_TRIVIAL(fatal) << "Fatal error! For a quick resolution keep a copy of input data!";

			}

			if (source_dataset != NULL) GDALClose(source_dataset);
			for (auto out_dataset = outlayers_datasets_map.begin(); out_dataset != outlayers_datasets_map.end(); ++out_dataset) {
				if (out_dataset->second != NULL) GDALClose(out_dataset->second);
			}
		}
	}
}