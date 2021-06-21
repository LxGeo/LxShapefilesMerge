#include "segment_graph.h"
#include "defs.h"
#include "defs_cgal.h"
#include "defs_boost.h"
#include "parameters.h"
#include "segment_graph_wrapper.h"
#include <assert.h>

namespace LxGeo
{
	namespace LxShapefilesMerge
	{
		SegmentGraph::SegmentGraph(std::vector<Inexact_Segment_2>& all_segments,
			std::vector<short int>& segment_LID,
			std::vector<short int>& segment_PID,
			std::vector<double>& segment_angle,
			Boost_RTree_2 segments_tree)
		{
			// load weights and thresh from parameters
			e_distance_weight = params->e_distance_weight;
			e_angle_weight = params->e_angle_weight;
			MAX_GROUPING_DISTANCE = params->MAX_GROUPING_DISTANCE;
			MAX_GROUPING_ANGLE_DEG = params->MAX_GROUPING_ANGLE_DEG;
			MAX_GROUPING_ANGLE_RAD = MAX_GROUPING_ANGLE_DEG * M_PI / 180;

			all_segments = all_segments;
			segment_LID = segment_LID;
			segment_PID = segment_PID;
			segment_angle = segment_angle;
			segments_tree = segments_tree;
			SG = BoostSegmentGraph(all_segments.size());
		}


		SegmentGraph::~SegmentGraph()
		{
		}

		void SegmentGraph::fill_graph()
		{

			for (size_t c_segment_index=0; c_segment_index <all_segments.size(); ++c_segment_index)
			{

				std::vector<size_t> neighborhood_indices;
				std::vector<double> neighborhood_weights;

				get_neighborhood_segments_indices_weights(c_segment_index,
					neighborhood_indices,
					neighborhood_weights
				);

				//test only
				assert (neighborhood_indices.size() == neighborhood_weights.size());

				double sum_vertex_centrality = 0;
				// adding graph edges
				for (size_t neigh_iter_index; neigh_iter_index < neighborhood_indices.size(); ++neigh_iter_index)
				{

					size_t neigh_index = neighborhood_indices[neigh_iter_index];
					double neigh_weight = neighborhood_weights[neigh_iter_index];
					
					auto e = boost::add_edge(c_segment_index, neigh_index, SG).first;

					// the weight can be assigned using a weight map.
					boost::property_map<BoostSegmentGraph, boost::edge_weight_t>::type weightmap =
						get(boost::edge_weight, SG);
					weightmap[e] = neigh_weight;
					sum_vertex_centrality += neigh_weight;

				}



				

			}

		}

		void SegmentGraph::get_neighborhood_segments_indices_weights(size_t c_segment_index,
			std::vector<size_t>& neighborhood_indices,
			std::vector<double>& neighborhood_weights
		)
		{
			Inexact_Segment_2 c_segment = all_segments[c_segment_index];
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

			segments_tree.query(bgi::intersects(query), std::back_inserter(possible_neighbors));

			for (int i = 0; i < possible_neighbors.size(); i++) {

				size_t possible_i_index = possible_neighbors[i].second;
				Inexact_Segment_2 possible_i_segment = all_segments[possible_i_index];

				Boost_Segment_2 possible_i_boost_segment = Boost_Segment_2(
					Boost_Point_2(possible_i_segment.source().x(), possible_i_segment.source().y()),
					Boost_Point_2(possible_i_segment.target().x(), possible_i_segment.target().y())
				);
				//check distance
				double diff_distance = bg::distance(c_boost_segment, possible_i_boost_segment);
				if (diff_distance < MAX_GROUPING_DISTANCE)
				{
					// check angle difference
					double diff_angle = segment_angle[possible_i_index] - segment_angle[c_segment_index];
					if (diff_angle < MAX_GROUPING_ANGLE_RAD)
					{
						//// ADD segment overlap ratios filter

						neighborhood_indices.push_back(possible_i_index);
						double edge_pre_weight = ((diff_distance / MAX_GROUPING_DISTANCE) + (diff_angle / MAX_GROUPING_ANGLE_RAD)) / 2;
						// below is to decrease edge weight based on LID & PID of segments
						if (segment_LID[c_segment_index] != segment_LID[possible_i_index]) neighborhood_weights.push_back(edge_pre_weight / 3);
						else if (segment_PID[c_segment_index] != segment_PID[possible_i_index]) neighborhood_weights.push_back(edge_pre_weight / 2);
						else neighborhood_weights.push_back(edge_pre_weight);
					}
				}
			}
		}
	}
}