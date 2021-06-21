#pragma once

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

struct VertexData
{

    size_t id;
    double centrality;
};


typedef boost::adjacency_list<boost::vecS, boost::vecS,
    boost::undirectedS,
    VertexData,
    boost::property<boost::edge_weight_t, double>
> BoostSegmentGraph;