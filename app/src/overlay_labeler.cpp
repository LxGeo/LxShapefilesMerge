#include "defs.h"
#include "defs_common.h"
#include "overlay_labeler.h"
#include "facet_graph_wrapper.h"
#include "boost/geometry/algorithms/point_on_surface.hpp"


namespace LxGeo
{
	namespace LxShapefilesMerge
	{
		using namespace GeometryFactoryShared;
		using namespace shapefiles_merge_utils;

		void OverlayLabeler::construct_graph() {

			//all_facets_polygons.reserve(all_facets.size());
			//auto facet_e_count_map = boost::get(&VertexData::facet_e_count, FG);
			//auto edge_e_layers_map = boost::get(&EdgeData::edge_e_layers, FG);

			// add background vertex with descriptor 0
			background_vertex = boost::add_vertex(FG);
			FG[background_vertex].facet_e_count = 0;
			// Iterate over all facets
			for (size_t facet_idx = 0; facet_idx < all_facets.size(); ++facet_idx) {

				// 
				const std::vector<std::vector<Inexact_Point_2>>& facet_polygon = all_facets[facet_idx];
				Boost_Polygon_2 b_facet_polygon;

				// construct boost polygon with holes
				for (auto c_ring_pt : facet_polygon[0]) bg::append(b_facet_polygon.outer(), transform_C2B_Point(c_ring_pt));
				/*for (size_t inner_ring_idx = 1; inner_ring_idx < facet_polygon.size(); ++inner_ring_idx) {
					for (auto c_ring_pt : facet_polygon[inner_ring_idx]) bg::append(b_facet_polygon.inners()[inner_ring_idx - 1], transform_C2B_Point(c_ring_pt));
				}*/
				// Get representative_point ( point that's within the polygon and outside its inner rings ) 
				
				double b_facet_area = bg::area(b_facet_polygon);
				if (b_facet_area < 1e-2) continue;
				
				Boost_Point_2 representative_point;
				boost::geometry::centroid(b_facet_polygon, representative_point); //point_on_surface
				size_t point_count = point_count_in_layers(representative_point);

				// add current facet to graph
				vertex_descriptor c_vertex = boost::add_vertex(FG);
				FG[c_vertex].facet_e_count = static_cast<short int>(point_count);

				// add boost_polygon to all_facets_polygons
				all_facets_polygons_map[c_vertex] = std::make_pair(b_facet_polygon, bg::area(b_facet_polygon));

				for (auto c_ring : facet_polygon) {

					for (int edge_idx = 0; edge_idx < c_ring.size() - 1; ++edge_idx) {
						// create current line_string
						Boost_LineString_2 c_edge_linestring;
						bg::append(c_edge_linestring, transform_C2B_Point(c_ring[edge_idx]));
						bg::append(c_edge_linestring, transform_C2B_Point(c_ring[edge_idx + 1]));

						vertex_descriptor source_vertex = c_vertex;
						vertex_descriptor target_vertex;
						std::vector<bool> layers_overlap_search_result;

						// search for spatially equal segment
						std::pair<bool, size_t> search_result = search_equal_linestring(c_edge_linestring);

						if (!search_result.first) {
							// if shared edge is new
							//// add shared edge
							add_to_shared_edges(c_edge_linestring);
							//// compute edge's layer overlap vector
							layers_overlap_search_result = get_shared_edge_layer_overlap(c_edge_linestring);
							//// Set target vertex to background
							target_vertex = background_vertex;

							auto add_pair = boost::add_edge(source_vertex, target_vertex, FG);
							FG[add_pair.first].edge_e_layers = layers_overlap_search_result;
							FG[add_pair.first].id = graph_edge_descriptors.size();
							FG[add_pair.first].length = bg::length(c_edge_linestring);
							graph_edge_descriptors.push_back(add_pair.first);
						}
						else {
							// if edge already added
							edge_descriptor c_edge = graph_edge_descriptors[search_result.second];
							layers_overlap_search_result = FG[c_edge].edge_e_layers;
							double edge_length = FG[c_edge].length;
							vertex_descriptor u = source(c_edge, FG);
							vertex_descriptor v = target(c_edge, FG);
							//if ((background_vertex != u) & (background_vertex != v)) throw std::exception("background_vertex check it !");

							target_vertex = (u == background_vertex) ? v : u;
							// remove existing edge
							boost::remove_edge(c_edge, FG);

							auto add_pair = boost::add_edge(source_vertex, target_vertex, FG);
							graph_edge_descriptors[search_result.second] =add_pair.first;
							FG[add_pair.first].edge_e_layers = layers_overlap_search_result;
							FG[add_pair.first].id = search_result.second;
							FG[add_pair.first].length = edge_length;
							if (!add_pair.second) throw std::exception("edge already added check it !");
						}

					}
				}


			}

			// convert all_shared_edges to all_ogr_shared_edges
			for (Boost_LineString_2 c_line : all_shared_edges) {
				OGRLineString converted_line;
				for (Boost_Point_2 c_point : c_line) {
					converted_line.addPoint(bg::get<0>(c_point), bg::get<1>(c_point));
				}
				all_ogr_shared_edges.push_back(converted_line);
			}
			//boost::clear_vertex(background_vertex, FG);
			boost::remove_vertex(background_vertex, FG);

		}

		size_t OverlayLabeler::point_count_in_layers(const Boost_Point_2& search_point) {

			size_t in = 0;
			Inexact_Point_2 cgal_search_point = transform_B2C_Point(search_point);
			size_t layers_count = polygon_trees.size();
			for (size_t c_layer_idx = 0; c_layer_idx < layers_count; ++c_layer_idx) {

				const Boost_RTree_2& RT = polygon_trees[c_layer_idx];

				Boost_LineString_2 temp_linestring;
				bg::append(temp_linestring, Boost_Point_2(bg::get<0>(search_point)+0.005, bg::get<1>(search_point) + 0.005));
				bg::append(temp_linestring, Boost_Point_2(bg::get<0>(search_point) - 0.005, bg::get<1>(search_point) - 0.005));
				Boost_Box_2 query;
				bg::envelope(temp_linestring, query);
				std::vector<Boost_Value_2> candidates;
				RT.query(bgi::intersects(query), std::back_inserter(candidates));

				bool includes = false;

				for (size_t u = 0; u < candidates.size(); ++u) {

					const std::vector<std::vector<Inexact_Point_2> >& P = all_regularized_polygons[c_layer_idx][candidates[u].second];

					// Tests inclusion of f_bar in P
					// Inclusion means : inside the outer ring of P, not inside the inner rings (holes) of P

					bool inside_outer_ring = false;
					bool inside_hole = false;
					Inexact_Polygon_2 Q_out(P[0].cbegin(), P[0].cend()-1);

					if (!Q_out.is_simple()) continue;
					inside_outer_ring = (Q_out.has_on_bounded_side(cgal_search_point));

					if (inside_outer_ring) {
						for (size_t v = 1; v < P.size(); ++v) {
							Inexact_Polygon_2 Q_in(P[v].cbegin(), P[v].cend());
							if (!Q_in.is_simple()) continue;
							inside_hole = (Q_in.has_on_bounded_side(cgal_search_point));
							if (inside_hole) break;
						}
					}

					if (inside_outer_ring && !inside_hole) {
						includes = true;
						break;
					}
				}

				// Final vote :
				if (includes) ++in;
			}

			return in;

		}

		std::vector<bool> OverlayLabeler::get_shared_edge_layer_overlap(const Boost_LineString_2& search_edge) {

			size_t layers_count = polygon_trees.size();
			std::vector<bool> search_result(layers_count, false);

			for (size_t c_layer_idx = 0; c_layer_idx < layers_count; ++c_layer_idx) {

				const std::vector<Boost_LineString_2>& c_layer_edges = all_regularized_edges[c_layer_idx];
				const Boost_RTree_2& c_layer_edges_trees = edges_trees[c_layer_idx];

				Boost_Box_2 query;
				bg::envelope(search_edge, query);
				std::vector<Boost_Value_2> candidates;
				c_layer_edges_trees.query(bgi::intersects(query), std::back_inserter(candidates));

				for (size_t u = 0; u < candidates.size(); ++u) {

					const Boost_LineString_2& c_candidate_linestring = c_layer_edges[candidates[u].second];
					Boost_LineString_2 intersection_output;
					boost::geometry::intersection(search_edge, c_candidate_linestring, intersection_output);
					//if (boost::geometry::covered_by(search_edge, c_candidate_linestring)) {
					if (intersection_output.size()>1){
						search_result[c_layer_idx] = true;
						break;
					}
				}
			}
			return search_result;
		}

		std::pair<bool, size_t> OverlayLabeler::search_equal_linestring(Boost_LineString_2 search_linestring) {

			std::pair<bool, size_t> search_result = { false, NULL };

			Boost_Box_2 query(search_linestring.at(0), search_linestring.at(1));
			std::vector<Boost_Value_2> candidates;
			shared_edges_tree.query(bgi::intersects(query), std::back_inserter(candidates));

			for (size_t u = 0; u < candidates.size(); ++u) {

				const Boost_LineString_2& c_candidate_linestring = all_shared_edges[candidates[u].second];
				if (boost::geometry::equals(search_linestring, c_candidate_linestring)) {
					search_result.first = true;
					search_result.second = candidates[u].second;
					break;
				}
			}

			return search_result;
		}

		void OverlayLabeler::add_to_shared_edges(const Boost_LineString_2& edge_to_add) {
			// append to all_shared_edges
			size_t new_edge_index = all_shared_edges.size();
			all_shared_edges.push_back(edge_to_add);

			Boost_Box_2 envelope;
			bg::envelope(edge_to_add, envelope);
			// add to all_shared_edges tree
			shared_edges_tree.insert(Boost_Value_2(envelope, new_edge_index));
		}

		std::vector<std::list<vertex_descriptor>> OverlayLabeler::connected_components_subgraphs()
		{
			/*vertex_component_map mapping = boost::make_shared<std::vector<unsigned long>>(num_vertices(FG));
			size_t num = boost::connected_components(FG, mapping->data());

			std::vector<ComponentGraph> component_graphs;

			for (size_t i = 0; i < num; i++)
				component_graphs.emplace_back(FG,
					[&](edge_descriptor e) {
						return mapping->at(boost::source(e, FG)) == i
							|| mapping->at(boost::target(e, FG)) == i;
					},
					[&](vertex_descriptor v) {
						return mapping->at(v) == i;
					});

			return component_graphs;*/

			vertex_descriptor_map idxMap;
			boost::associative_property_map<vertex_descriptor_map> indexMap(idxMap);

			vertex_iterator di, dj;
			boost::tie(di, dj) = boost::vertices(FG);
			for (int i = 0; di != dj; ++di, ++i) {
				boost::put(indexMap, (*di), i);
			}

			std::map<vertex_descriptor, size_t> compMap;
			boost::associative_property_map<vertex_descriptor_map> componentMap(compMap);

			size_t num = boost::connected_components(FG, componentMap, boost::vertex_index_map(indexMap));

			std::vector<std::list<vertex_descriptor>> component_graphs(num);
			for (auto comp_element : compMap) {
				component_graphs[comp_element.second].push_back(comp_element.first);
			}

			
			return component_graphs;
		}

		void OverlayLabeler::log_subgraphs() {

			std::vector<std::list<vertex_descriptor>> component_sub_graphs = connected_components_subgraphs();
			for (auto const& component : component_sub_graphs) {		
				extract_geometries_from_subgraph3(component);				
			}

			for (auto poly_term_pair : all_extracted_polygons_terms) {
				std::cout << "Score: " << poly_term_pair.second << std::endl;
				std::cout << poly_term_pair.first.exportToWkt() << std::endl;
			}
		}

		void OverlayLabeler::get_edge_partition(const std::set<vertex_descriptor>& selected_facets,
			std::set<edge_descriptor>& used_edges, std::set<edge_descriptor>& ignored_edges,
			std::set<vertex_descriptor>& current_selection_neighbours) {

			// clear used & ignored edges
			used_edges.clear();
			ignored_edges.clear();
			current_selection_neighbours.clear();

			for (vertex_descriptor c_selected_facet : selected_facets) {
				out_edge_iterator e, e_end;
				for (boost::tie(e, e_end) = boost::out_edges(c_selected_facet, FG); e != e_end; ++e) {

					vertex_descriptor edge_source_v = boost::source(*e, FG);
					bool source_in_selected = (selected_facets.find(edge_source_v) != selected_facets.end());
					if (!source_in_selected) current_selection_neighbours.insert(edge_source_v);
					vertex_descriptor edge_target_v = boost::target(*e, FG);
					bool target_in_selected = (selected_facets.find(edge_target_v) != selected_facets.end());
					if (!target_in_selected) current_selection_neighbours.insert(edge_target_v);
					
					if (source_in_selected && target_in_selected) {
						ignored_edges.insert(*e);
					}
					else {
						used_edges.insert(*e);
					}

				}
			}
		}

		double OverlayLabeler::compute_edge_diff_term(const std::set<edge_descriptor>& used_edges, const std::set<edge_descriptor>& ignored_edges) {

			double all_edges_length = 0;
			double opt_term = 0;

			// used edges iteration
			for (const edge_descriptor& c_edge_desc : used_edges) {
				auto c_edge_map = FG[c_edge_desc];
				double c_edge_length = c_edge_map.length; // this should be changed (length must be computed once at graph construction)
				std::vector<bool> c_edge_e_layers = c_edge_map.edge_e_layers;
				std::vector<int> c_edge_e_layers_value(c_edge_e_layers.begin(), c_edge_e_layers.end());

				double c_edge_term = (inner_product(
					c_edge_e_layers_value.begin(),
					c_edge_e_layers_value.end(),
					edge_diff_layer_weights.begin(), double(0)) / double(edge_diff_layer_weights_sum));// -0.5;
				opt_term = opt_term + c_edge_term * c_edge_length;
				all_edges_length += c_edge_length;
			}

			// ignored edges iteration
			for (const edge_descriptor& c_edge_desc : ignored_edges) {
				auto c_edge_map = FG[c_edge_desc];
				double c_edge_length = c_edge_map.length;
				std::vector<bool> c_edge_e_layers = c_edge_map.edge_e_layers;
				std::vector<int> c_edge_e_layers_value(c_edge_e_layers.begin(), c_edge_e_layers.end());

				double c_edge_term = (inner_product(
					c_edge_e_layers_value.begin(),
					c_edge_e_layers_value.end(),
					edge_diff_layer_weights.begin(), double(0)) / double(edge_diff_layer_weights_sum));// -0.5;
				opt_term = opt_term - c_edge_term * c_edge_length;
				all_edges_length += c_edge_length;
			}

			//normalize bu totla length
			opt_term /= all_edges_length;
			return opt_term;

		}

		void OverlayLabeler::reconstruct_polygon_from_edges(const std::set<edge_descriptor>& used_edges, std::string& polygon_to_reconstruct) {

			std::vector<Boost_LineString_2> used_linestrings(used_edges.size());
			std::transform(used_edges.begin(), used_edges.end(), used_linestrings.begin(),
				[this](edge_descriptor c_ed_desc)->Boost_LineString_2 { return all_shared_edges[FG[c_ed_desc].id]; }
			);

			OGRMultiLineString linestring_work;
			for (const edge_descriptor& c_edge_desc : used_edges) {
				const OGRGeometry& c_ogr_linestring= all_ogr_shared_edges[FG[c_edge_desc].id];
				linestring_work.addGeometry(&c_ogr_linestring);
			}

			OGRGeometry *polygonized= linestring_work.Polygonize();
			std::cout << polygonized->exportToWkt()<< std::endl;
			polygon_to_reconstruct = polygonized->exportToWkt();
			//polygon_to_reconstruct = *polygonized->toPolygon();

		}
		
		/*void OverlayLabeler::extract_geometries_from_subgraph(std::vector<vertex_descriptor> subgraph_vertices) {

			// step0: transform subgraph vertices to map with keys as vertex descriptor and values a pair of ( size_t: facet_e_count && double: facet area)
			std::map<vertex_descriptor, std::pair<size_t, double>> available_facets_map;
			size_t max_facet_e_count = 1;
			// step0: define extracted polygons container
			std::list<OGRPolygon> extracted_polygons;
			std::list<double> extracted_polygons_opt_value;
			
			// step1: Get all facets facet_e_count
			for (vertex_descriptor v_d : subgraph_vertices) {
				size_t c_facet_e_count = FG[v_d].facet_e_count;
				if (c_facet_e_count > max_facet_e_count) max_facet_e_count = c_facet_e_count;
				double c_facet_area = bg::area(all_facets_polygons[size_t(v_d)]);
				available_facets_map[v_d] = std::pair(c_facet_e_count, c_facet_area);
			}

			// step2: initialize optimization vars
			double c_selection_opt_value = -1;
			std::set<vertex_descriptor> c_selected_facets;
			OGRPolygon last_polygon_reconstructed;
			// step3: selection loop
			while (!available_facets_map.empty()) {
				//step a: push first facet c_selected_facets 
				vertex_descriptor last_added_vertex = available_facets_map.begin()->first;
				c_selected_facets.insert(last_added_vertex);

				std::set<edge_descriptor> used_edges, ignored_edges;
				get_edge_partition(c_selected_facets, used_edges, ignored_edges);

				OGRPolygon selection_reconstructed_polygon;
				reconstruct_polygon_from_edges(used_edges, selection_reconstructed_polygon);

				double edges_diff_term = compute_edge_diff_term(used_edges, ignored_edges);

				if (edges_diff_term < c_selection_opt_value) {
					c_selection_opt_value = 0;
					extracted_polygons.push_back(last_polygon_reconstructed);
					extracted_polygons_opt_value.push_back(edges_diff_term);
					//for (vertex_descriptor c_chosen_vertex_desc : c_selected_facets) { available_facets_map.erase(c_chosen_vertex_desc); }
					c_selected_facets.clear();
				}
				else {
					last_polygon_reconstructed = selection_reconstructed_polygon;
					c_selection_opt_value = edges_diff_term;
					available_facets_map.erase(last_added_vertex);
				}

			}
			if (!c_selected_facets.empty()) {
				extracted_polygons.push_back(last_polygon_reconstructed);
				extracted_polygons_opt_value.push_back(c_selection_opt_value);
			}

			all_extracted_polygon.insert(all_extracted_polygon.end(), extracted_polygons.begin(), extracted_polygons.end());
			all_extracted_polygon_opt_values.insert(all_extracted_polygon_opt_values.end(), extracted_polygons_opt_value.begin(), extracted_polygons_opt_value.end());
			

		}*/

		/*void OverlayLabeler::extract_geometries_from_subgraph2(std::list<vertex_descriptor> subgraph_vertices) {

			std::list<OGRPolygon> extracted_polygons;
			std::list<double> extracted_polygons_opt_term;
			std::map<vertex_descriptor, std::pair<size_t, double>> available_facets_map;
			
			for (vertex_descriptor v_d : subgraph_vertices) {
				size_t c_facet_e_count = FG[v_d].facet_e_count;
				double c_facet_area = all_facets_polygons_map[v_d].second;
				available_facets_map[v_d] = std::pair(c_facet_e_count, c_facet_area);
			}
			// create a vector of vertices from map (ordered by facet_e_count)
			std::list<vertex_descriptor> orderd_vertices_descriptors;
			for (auto it = available_facets_map.begin(); it != available_facets_map.end(); ++it) {orderd_vertices_descriptors.push_back(it->first);}
			orderd_vertices_descriptors.sort(
				[&](const vertex_descriptor& a, const vertex_descriptor& b) -> bool {return available_facets_map[a].first > available_facets_map[b].first; }
			);


			// sets describing current selection edge partition
			std::set<edge_descriptor> used_edges, ignored_edges;
			std::pair<vertex_descriptor, std::pair<size_t, double>> last_chosen_KV;
			// get first in ordered (facet with highest e_count)
			vertex_descriptor starting_vertex = orderd_vertices_descriptors.front();
			std::set<vertex_descriptor> c_selection_facets = { starting_vertex };
			// each confirmed addition to selection must followed by item removal from both map and ordered_vertices
			available_facets_map.erase(starting_vertex);
			orderd_vertices_descriptors.pop_front();
			//last_chosen_KV = std::make_pair( starting_vertex, available_facets_map[starting_vertex] );

			std::set<vertex_descriptor> current_selection_neighbours;
			//update current_selection state variables
			get_edge_partition(c_selection_facets, used_edges, ignored_edges, current_selection_neighbours);
			// remove neighbours that are not available
			for (auto c_neighbour = current_selection_neighbours.begin(); c_neighbour != current_selection_neighbours.end(); ) {
				if (available_facets_map.find(*c_neighbour) == available_facets_map.end())
					c_neighbour = current_selection_neighbours.erase(c_neighbour);
				else ++c_neighbour;
			}
			
			double last_opt_term = 0;
			OGRPolygon last_selection_polygon;
			while (!current_selection_neighbours.empty()) {

				double edges_diff_term = compute_edge_diff_term(used_edges, ignored_edges);

				// augment selection
				if (last_opt_term < edges_diff_term) {
					// choose next facet to add from neighbours ( neighbour with lowest facet_e_count )
					vertex_descriptor chosen_neighbour;
					for (auto c_possible_neighbour = orderd_vertices_descriptors.rbegin(); c_possible_neighbour != orderd_vertices_descriptors.rend(); ++c_possible_neighbour) {
						if (current_selection_neighbours.find(*c_possible_neighbour) != current_selection_neighbours.end()) {
							chosen_neighbour = *c_possible_neighbour;
							break;
						}
					}

					// updat last selection variables
					last_opt_term = edges_diff_term;
					reconstruct_polygon_from_edges(used_edges, last_selection_polygon);

					last_chosen_KV = { chosen_neighbour, available_facets_map[chosen_neighbour] };
					c_selection_facets.insert(chosen_neighbour);
					available_facets_map.erase(starting_vertex);
					orderd_vertices_descriptors.remove(chosen_neighbour);
					get_edge_partition(c_selection_facets, used_edges, ignored_edges, current_selection_neighbours);
					// remove neighbours that are not available
					for (auto c_neighbour = current_selection_neighbours.begin(); c_neighbour != current_selection_neighbours.end(); ) {
						if (available_facets_map.find(*c_neighbour) == available_facets_map.end())
							c_neighbour = current_selection_neighbours.erase(c_neighbour);
						else ++c_neighbour;
					}
				}
				// readd last removed facet to avaialbale // add last polygon with opt_term // init new selection // init current_selection_neighbours
				else {
					available_facets_map[last_chosen_KV.first] = last_chosen_KV.second;
					orderd_vertices_descriptors.push_back(last_chosen_KV.first);
					// resort vertices_descriptors list after insertion
					orderd_vertices_descriptors.sort(
						[&](const vertex_descriptor& a, const vertex_descriptor& b) -> bool {return available_facets_map[a].first > available_facets_map[b].first; }
					);

					if (current_selection_neighbours.empty()){
					// this case is important
					}
					extracted_polygons.push_back(last_selection_polygon);
					last_selection_polygon = OGRPolygon();
					extracted_polygons_opt_term.push_back(last_opt_term);
					last_opt_term = 0;

					// create new selection
					vertex_descriptor starting_vertex = orderd_vertices_descriptors.front();
					c_selection_facets = { starting_vertex };
					//update current_selection state variables
					get_edge_partition(c_selection_facets, used_edges, ignored_edges, current_selection_neighbours);
					// remove neighbours that are not available
					for (auto  c_neighbour = current_selection_neighbours.begin(); c_neighbour != current_selection_neighbours.end(); ) {
						if (available_facets_map.find(*c_neighbour) == available_facets_map.end())
							c_neighbour = current_selection_neighbours.erase(c_neighbour);
						else ++c_neighbour;
					}
				}

			}

			// add last polygon
			if (last_opt_term > 0) { //must be changed
				double edges_diff_term = compute_edge_diff_term(used_edges, ignored_edges);
				last_opt_term = edges_diff_term;

				reconstruct_polygon_from_edges(used_edges, last_selection_polygon);
				extracted_polygons.push_back(last_selection_polygon);
				extracted_polygons_opt_term.push_back(last_opt_term);
			}

			all_extracted_polygon.insert(all_extracted_polygon.end(), extracted_polygons.begin(), extracted_polygons.end());
			all_extracted_polygon_opt_values.insert(all_extracted_polygon_opt_values.end(), extracted_polygons_opt_term.begin(), extracted_polygons_opt_term.end());
		}*/

		void OverlayLabeler::extract_geometries_from_subgraph3(std::list<vertex_descriptor> subgraph_vertices) {

			std::list<std::pair<std::string, double>> extracted_polygons_terms;
			std::map<vertex_descriptor, std::pair<size_t, double>> available_facets_map;

			for (vertex_descriptor v_d : subgraph_vertices) {
				if (v_d == background_vertex) continue;
				size_t c_facet_e_count = FG[v_d].facet_e_count;
				double c_facet_area = all_facets_polygons_map[v_d].second;
				available_facets_map[v_d] = std::pair(c_facet_e_count, c_facet_area);
			}
			if (available_facets_map.size() == 0)return;

			// state variable
			std::set<vertex_descriptor> c_selection_set;
			double last_opt_term;
			std::string last_polygon;
			//std::queue<vertex_descriptor> possible_neighbours_queue;
			std::set<vertex_descriptor> possible_neighbours_set;


			auto get_vertex_neighbours = [&](const vertex_descriptor& c_vertex, const std::set<vertex_descriptor>& ignore_set) {

				std::set<vertex_descriptor> neighbours_set;
				out_edge_iterator e, e_end;
				for (boost::tie(e, e_end) = boost::out_edges(c_vertex, FG); e != e_end; ++e) {
					vertex_descriptor edge_source_v = boost::source(*e, FG);
					bool source_in_selected = (ignore_set.find(edge_source_v) != ignore_set.end());
					bool source_in_available = (available_facets_map.find(edge_source_v) != available_facets_map.end());
					if (!source_in_selected & source_in_available) neighbours_set.insert(edge_source_v);
					vertex_descriptor edge_target_v = boost::target(*e, FG);
					bool target_in_selected = (ignore_set.find(edge_target_v) != ignore_set.end());
					bool target_in_available = (available_facets_map.find(edge_target_v) != available_facets_map.end());
					if (!target_in_selected & target_in_available) neighbours_set.insert(edge_target_v);
				}
				if (neighbours_set.find(background_vertex) != neighbours_set.end()) neighbours_set.erase(background_vertex);
				return neighbours_set;
			};

			// lambda to add to current selection & remove from available & adds respective to possible_neighbours
			auto add_to_current_selection = [&](vertex_descriptor vertex_to_add) {
				c_selection_set.insert(vertex_to_add);
				available_facets_map.erase(vertex_to_add);

				std::set<vertex_descriptor> list_of_neighbours = get_vertex_neighbours(vertex_to_add, c_selection_set);
				for (vertex_descriptor c_poss_neigh : list_of_neighbours) possible_neighbours_set.insert(c_poss_neigh);
			};

			// lambda to get facet with heighest or lowest facet_e_count from available_facets_map
			auto get_low_high_facet = [&available_facets_map](bool heighest)
			{
				auto result = std::max_element(std::begin(available_facets_map), std::end(available_facets_map),
					[&heighest](const auto a, const auto b)
					{
						if (heighest)
							return a.second.first < b.second.first;
						else
							return a.second.first > b.second.first;
					});
				return result->first;
			};

			auto construct_from_selection = [&](const std::set<vertex_descriptor> selection_set) {
				std::set<edge_descriptor> used_edges, ignored_edges;
				std::set<vertex_descriptor> neighbours_set; // should be ignored
				get_edge_partition(selection_set, used_edges, ignored_edges, neighbours_set);
				double edges_diff_term = compute_edge_diff_term(used_edges, ignored_edges);
				std::string constructed_polygon;
				std::cout << "term: " << edges_diff_term << std::endl;
				reconstruct_polygon_from_edges(used_edges, constructed_polygon);
				return std::make_pair( constructed_polygon, edges_diff_term );
			};

			vertex_descriptor start_vertex = get_low_high_facet(true);
			add_to_current_selection(start_vertex);
			std::pair<std::string, double> polygon_term_pair = construct_from_selection(c_selection_set);
			last_polygon = polygon_term_pair.first; last_opt_term = polygon_term_pair.second;

			while (!c_selection_set.empty()) {

				if (possible_neighbours_set.empty()) {
					std::cout << "starting new: " << std::endl;
					// push last created polygon
					extracted_polygons_terms.push_back(std::make_pair(last_polygon, last_opt_term));
					// clear state_variable
					c_selection_set.clear();
					// if still facets create new starting vertex
					if (!available_facets_map.empty()) {
						//start new selection
						vertex_descriptor start_vertex = get_low_high_facet(true);
						add_to_current_selection(start_vertex);
						std::pair<std::string, double> polygon_term_pair = construct_from_selection(c_selection_set);
						last_polygon = polygon_term_pair.first; last_opt_term = polygon_term_pair.second;
					}
				}
				else {
					// pop from set == get first and erase
					vertex_descriptor c_possible_neighbour = *possible_neighbours_set.begin();
					possible_neighbours_set.erase(c_possible_neighbour);
					// temporary selection using c_possible neighbour
					std::set<vertex_descriptor> temporary_selection(c_selection_set.begin(), c_selection_set.end());
					temporary_selection.insert(c_possible_neighbour);

					std::pair<std::string, double> polygon_term_pair = construct_from_selection(temporary_selection);
					if (polygon_term_pair.second > last_opt_term) {
						// update last state
						last_opt_term = polygon_term_pair.second; last_polygon = polygon_term_pair.first;
						add_to_current_selection(c_possible_neighbour);
					}
					else std::cout<<"Ignored"<<std::endl;

				}

			}
		
			//addition to global polygons
			//all_extracted_polygons_terms.insert(all_extracted_polygons_terms.end(), extracted_polygons_terms.begin(), extracted_polygons_terms.end());

		}

	}

}