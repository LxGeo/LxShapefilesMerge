#pragma once
#include "defs.h"
#include "defs_cgal.h"
#include "defs_boost.h"
#include "segment_graph_wrapper.h"

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
			SegmentGraph(std::vector<Inexact_Segment_2>& all_segments,
				std::vector<short int>& segment_LID,
				std::vector<short int>& segment_PID,
				std::vector<double>& segment_angle,
				Boost_RTree_2& segments_tree);

			~SegmentGraph();

			void fill_graph();

			void SegmentGraph::get_neighborhood_segments_indices_weights(size_t c_segment_index,
				std::vector<size_t>& neighborhood_indices,
				std::vector<double>& neighborhood_weights,
				std::vector<double>& neighborhood_distance
			);

			void SegmentGraph::cluster_segments();

			size_t SegmentGraph::get_min_centrality_vertex_index();

			void SegmentGraph::get_connected_vertices_indices(size_t vertex_idx, std::vector<size_t>& c_connected_vertices_indices);

			void SegmentGraph::create_segment_group(size_t c_group_id, std::vector<size_t> c_connected_vertices_indices);

			void SegmentGraph::delete_vertex_related(size_t vertex_idx);

		public:
			BoostSegmentGraph SG;
		private:
			double e_distance_weight;
			double e_angle_weight;
			double MAX_GROUPING_DISTANCE;
			double MAX_GROUPING_ANGLE_DEG;
			double MAX_GROUPING_ANGLE_RAD;
			std::vector<Inexact_Segment_2> _all_segments;
			std::vector<short int>& _segment_LID;
			std::vector<short int>& _segment_PID;
			std::vector<double> _segment_angle;
			Boost_RTree_2 _segments_tree;
			size_t groupes_count;
			std::vector<double> vertcies_centrality;
			std::vector<size_t> vertcies_groups;
		};
	}
}