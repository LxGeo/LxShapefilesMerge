#pragma once
#include "defs.h"
#include "defs_cgal.h"
#include "defs_boost.h"
#include "segment_graph_wrapper.h"
#include <ogrsf_frmts.h>

namespace LxGeo
{
	namespace LxShapefilesMerge
	{
		using namespace GeometryFactoryShared;
		/**
		*  A SegmentGraph class to group (cluster) segments.
		*/
		class SegmentGraph
		{
		public:
			SegmentGraph() = default;
			SegmentGraph(std::vector<Segment_2>& all_segments,
				std::vector<short int>& segment_LID,
				std::vector<short int>& segment_PID,
				std::vector<short int>& segment_ORDinP,
				std::vector<short int>& segment_RRSize,
				std::vector<double>& segment_angle,
				Boost_RTree_2& segments_tree);

			~SegmentGraph();

			void fill_graph();

			void SegmentGraph::get_neighborhood_segments_indices_weights(size_t c_segment_index,
				std::vector<size_t>& neighborhood_indices,
				std::vector<double>& neighborhood_weights,
				std::vector<double>& neighborhood_distance
			);

			void SegmentGraph::fix_n2_consecutive_segments();

			void SegmentGraph::cluster_segments();

			size_t SegmentGraph::get_min_centrality_vertex_index();

			void SegmentGraph::get_connected_vertices_indices(size_t vertex_idx, std::vector<size_t>& c_connected_vertices_indices);

			void SegmentGraph::filter_same_polygon_adjacent_indices(std::list<size_t>& c_connected_vertices_indices);

			size_t SegmentGraph::create_segment_group(size_t c_group_id, std::vector<size_t> c_connected_vertices_indices);

			void SegmentGraph::delete_vertex_related(size_t vertex_idx);

			void SegmentGraph::write_grouped_segments_shapefile( const std::string& output_filename);

			void SegmentGraph::fuse_segments();

			void SegmentGraph::get_respective_segment(size_t groupe_idx, std::vector<Segment_2*>& respective_segments);

			void SegmentGraph::get_best_fitting_line(Line_2& fitted_line, std::vector<Segment_2*>& respective_segments);
			
			void SegmentGraph::get_best_fitting_line(Line_2& fitted_line, std::vector<Segment_2*>& respective_segments, std::vector<size_t>& segments_weights);

			void SegmentGraph::get_best_fitting_line_by_direction(Line_2& fitted_line, std::vector<Segment_2*>& respective_segments);

			void SegmentGraph::get_best_fitting_line_by_direction(Line_2& fitted_line, std::vector<Segment_2*>& respective_segments, std::vector<size_t>& segments_weights);

			void SegmentGraph::add_polygon_to_layer(std::list<OGRLinearRing>& ex_int_rings, OGRLayer* current_dataset_layer, size_t current_segment_index);

			void SegmentGraph::reconstruct_polygons(const std::string& temp_dir);

		public:
			BoostSegmentGraph SG = BoostSegmentGraph(10);
		private:
			double e_distance_weight;
			double e_angle_weight;
			double MAX_GROUPING_DISTANCE;
			double MAX_GROUPING_ANGLE_DEG;
			double MAX_GROUPING_ANGLE_RAD;
			std::vector<Segment_2> _all_segments;
			std::vector<short int>& _segment_LID;
			std::vector<short int>& _segment_PID;
			std::vector<short int>& _segment_ORDinP;
			std::vector<short int>& _segment_RRSize;
			std::vector<double> _segment_angle;
			Boost_RTree_2 _segments_tree;
			size_t groupes_count;
			std::vector<double> vertcies_centrality;
			std::map<size_t, std::vector<size_t>> groups_map;
		public:
			std::vector<size_t> vertcies_groups;
		};
	}
}