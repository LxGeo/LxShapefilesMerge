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

			//auto facet_e_count_map = boost::get(&VertexData::facet_e_count, FG);
			//auto edge_e_layers_map = boost::get(&EdgeData::edge_e_layers, FG);

			// add background vertex with descriptor 0
			vertex_descriptor background_vertex = boost::add_vertex(FG);
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
				Boost_Point_2 representative_point;
				boost::geometry::centroid(b_facet_polygon, representative_point); //point_on_surface
				size_t point_count = point_count_in_layers(representative_point);

				// add current facet to graph
				vertex_descriptor c_vertex = boost::add_vertex(FG);
				FG[c_vertex].facet_e_count = static_cast<short int>(point_count);

				for (auto c_ring : facet_polygon) {

					for (size_t edge_idx = 0; edge_idx < c_ring.size() - 1; ++edge_idx) {
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
							graph_edge_descriptors.push_back(add_pair.first);
						}
						else {
							// if edge already added
							edge_descriptor c_edge = graph_edge_descriptors[search_result.second];
							layers_overlap_search_result = FG[c_edge].edge_e_layers;
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
							if (!add_pair.second) throw std::exception("edge already added check it !");
						}

					}
				}


			}

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
					Inexact_Polygon_2 Q_out(P[0].cbegin(), P[0].cend());

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
					if (intersection_output.size()>0){
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
	}

}