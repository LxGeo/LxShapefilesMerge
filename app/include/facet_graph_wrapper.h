#pragma once
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_matrix.hpp>

struct VertexData
{
    short int facet_e_count; // value represents count of facet in different layers
};

struct EdgeData
{
    std::vector<bool> edge_e_layers; // vector of edge existence in layers
};

typedef boost::adjacency_list<boost::vecS, boost::vecS,
    boost::undirectedS,
    boost::property<boost::vertex_index_t, size_t, VertexData>,
    boost::property<boost::edge_index_t, size_t, EdgeData>
> BoostFacetGraph;

//typedef boost::adjacency_matrix<boost::undirectedS,
//    boost::property<boost::vertex_index_t, size_t, VertexData>,
//    boost::property<boost::edge_index_t, size_t, EdgeData>
//>BoostFacetGraph;

typedef boost::graph_traits<BoostFacetGraph>::vertex_iterator vertex_iterator;
typedef boost::graph_traits<BoostFacetGraph>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<BoostFacetGraph>::edge_descriptor edge_descriptor;
typedef boost::graph_traits<BoostFacetGraph>::adjacency_iterator adjacency_iterator;
