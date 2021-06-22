#pragma once

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

struct VertexData
{
    double centrality;
};

struct EdgeData
{
    double distance;
};


typedef boost::adjacency_list<boost::vecS, boost::vecS,
    boost::undirectedS,
    //VertexData,
    boost::property<boost::vertex_index_t, size_t, boost::property< boost::vertex_centrality_t, double> >,
    boost::property<boost::edge_weight_t, double, EdgeData>
> BoostSegmentGraph;

typedef boost::graph_traits<BoostSegmentGraph>::vertex_iterator vertex_iterator;
typedef boost::graph_traits<BoostSegmentGraph>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<BoostSegmentGraph>::adjacency_iterator adjacency_iterator;