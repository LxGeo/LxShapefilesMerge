#pragma once
#include "defs.h"
#include "defs_boost.h"
#include "defs_cgal.h"
#include "shapefiles_merge_utilities.h"
#include "facet_graph_wrapper.h"
#include "parameters.h"
#include <thread>
#include <mutex>

namespace LxGeo
{
	namespace LxShapefilesMerge
	{
		using namespace GeometryFactoryShared;
		using namespace shapefiles_merge_utils;

		class OverlayLabeler
		{
		public:
			OverlayLabeler() {};

			OverlayLabeler(std::string overlay_layer_path,
				std::vector<std::string> regularized_layers_paths) {

				//load regularizd polygons
				read_shapefiles(regularized_layers_paths, all_regularized_polygons);
				// Extract edges from regularized polygons
				all_regularized_edges.resize(all_regularized_polygons.size());
				for (size_t layer_idx = 0; layer_idx < all_regularized_polygons.size(); ++layer_idx) {
					extract_edges_from_polygons(all_regularized_polygons[layer_idx], all_regularized_edges[layer_idx]);
				}

				//// generate rTrees for each layer of regularized (polygons and edges)
				polygon_trees.resize(all_regularized_polygons.size());
				edges_trees.resize(all_regularized_polygons.size());
				for (size_t layer_idx = 0; layer_idx < regularized_layers_paths.size(); ++layer_idx) {
					make_rtree_polygons(all_regularized_polygons[layer_idx], polygon_trees[layer_idx]);
					make_rtree_linestrings(all_regularized_edges[layer_idx], edges_trees[layer_idx]);
				}

				//load overlay facets
				read_single_shapefile(overlay_layer_path, NULL, all_facets);

				// load terms weights
				edge_diff_layer_weights = params->edge_diff_layer_weights;
				edge_diff_layer_weights_sum = std::accumulate(edge_diff_layer_weights.begin(), edge_diff_layer_weights.end(), size_t(0));

				//FG = BoostFacetGraph(all_facets.size()+1);
			};
			
			void construct_graph();

			~OverlayLabeler() {};

			size_t point_count_in_layers(const Boost_Point_2& search_point);

			std::vector<bool> get_shared_edge_layer_overlap(const Boost_LineString_2& search_edge);

			std::pair<bool, size_t> search_equal_linestring(Boost_LineString_2 search_linestring);

			void add_to_shared_edges(const Boost_LineString_2& search_edge);

			std::vector<std::list<vertex_descriptor>> connected_components_subgraphs();

			void extract_geometries_from_subgraph(std::vector<vertex_descriptor> subgraph_vertices);

			void log_subgraphs();

			void get_edge_partition(const std::set<vertex_descriptor>& selected_facets,
				std::set<edge_descriptor>& used_edges, std::set<edge_descriptor>& ignored_edges,
				std::set<vertex_descriptor>& current_selection_neighbours);

			void reconstruct_polygon_from_edges(const std::set<edge_descriptor>& used_edges, std::string& polygon_to_reconstruct);

			double compute_edge_diff_term(const std::set<edge_descriptor>& used_edges, const std::set<edge_descriptor>& ignored_edges);

			void extract_geometries_from_subgraph2(std::list<vertex_descriptor> subgraph_vertices);

			void extract_geometries_from_subgraph3(std::list<vertex_descriptor> subgraph_vertices);



		public:
			BoostFacetGraph FG;
			vertex_descriptor background_vertex;
			std::map<vertex_descriptor, int> subclusters;
			//std::vector<OGRPolygon> all_extracted_polygon;
			//std::vector<double> all_extracted_polygon_opt_values;
			std::list<std::pair<OGRPolygon, double>> all_extracted_polygons_terms;
		// overlay layer related attributes
		private:
			std::vector<std::vector<std::vector<Inexact_Point_2>>> all_facets; // convert later to cgal to use edge_iterator && get vertcies type
			std::vector<Boost_LineString_2> all_shared_edges; //used boost linestring to use equals (spatial) and overlaps
			std::vector<OGRLineString> all_ogr_shared_edges;
			Boost_RTree_2 shared_edges_tree;
			std::map<vertex_descriptor, std::pair<Boost_Polygon_2, double>> all_facets_polygons_map;//std::vector<Boost_Polygon_2> all_facets_polygons;
		//regularized layers related attributes
		private:
			std::vector<std::vector<std::vector<std::vector<Inexact_Point_2>>>> all_regularized_polygons; // 
			std::vector<Boost_RTree_2> polygon_trees; // regularized polygons Rtrees
			std::vector<std::vector<Boost_LineString_2>> all_regularized_edges;
			std::vector<Boost_RTree_2> edges_trees; // regularized edges Rtrees

			boost::property_map<BoostFacetGraph, boost::vertex_index_t>::type vertex_idx_map = boost::get(boost::vertex_index, FG);
			boost::property_map<BoostFacetGraph, boost::edge_index_t>::type edge_idx_map = boost::get(boost::edge_index, FG);
			std::vector<edge_descriptor> graph_edge_descriptors;

		private:
			std::mutex m_mutex;
			std::vector<size_t> edge_diff_layer_weights;
			size_t edge_diff_layer_weights_sum;

		};
	}
}