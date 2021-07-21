#pragma once
#include "defs.h"
#include "defs_boost.h"
#include "defs_cgal.h"
#include "shapefiles_merge_utilities.h"
#include "facet_graph_wrapper.h"

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
				all_regularized_edges.reserve(all_regularized_polygons.size());
				for (size_t layer_idx = 0; layer_idx < all_regularized_polygons.size(); ++layer_idx) {
					extract_edges_from_polygons(all_regularized_polygons[layer_idx], all_regularized_edges[layer_idx]);
				}

				//// generate rTrees for each layer of regularized (polygons and edges)
				for (size_t layer_idx = 0; layer_idx < regularized_layers_paths.size(); ++layer_idx) {
					make_rtree_polygons(all_regularized_polygons[layer_idx], polygon_trees[layer_idx]);
					make_rtree_linestrings(all_regularized_edges[layer_idx], edges_trees[layer_idx]);
				}

				//load overlay facets
				read_single_shapefile(overlay_layer_path, NULL, all_facets);

				FG = BoostFacetGraph(all_facets.size()+1);
			};
			
			void OverlayLabeler::construct_graph();

			~OverlayLabeler() {};

			size_t OverlayLabeler::point_count_in_layers(const Boost_Point_2& search_point);

			std::vector<bool> OverlayLabeler::get_shared_edge_layer_overlap(const Boost_LineString_2& search_edge);

			std::pair<bool, size_t> OverlayLabeler::search_equal_linestring(Boost_LineString_2 search_linestring);

			void OverlayLabeler::add_to_shared_edges(const Boost_LineString_2& search_edge);


		public:
			BoostFacetGraph FG;
		// overlay layer related attributes
		private:
			std::vector<std::vector<std::vector<Inexact_Point_2>>> all_facets; // convert later to cgal to use edge_iterator && get vertcies type
			std::vector<Boost_LineString_2> all_shared_edges; //used boost linestring to use equals (spatial) and overlaps
			Boost_RTree_2 shared_edges_tree;
		//regularized layers related attributes
		private:
			std::vector<std::vector<std::vector<std::vector<Inexact_Point_2>>>> all_regularized_polygons; // 
			std::vector<Boost_RTree_2> polygon_trees; // regularized polygons Rtrees
			std::vector<std::vector<Boost_LineString_2>> all_regularized_edges;
			std::vector<Boost_RTree_2> edges_trees; // regularized edges Rtrees

			boost::property_map<BoostFacetGraph, boost::vertex_index_t>::type vertex_idx_map = boost::get(boost::vertex_index_t(), FG);
			boost::property_map<BoostFacetGraph, boost::edge_index_t>::type edge_idx_map = boost::get(boost::edge_index_t(), FG);

		};
	}
}