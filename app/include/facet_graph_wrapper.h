#pragma once
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_matrix.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/make_shared.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/function.hpp>

struct VertexData
{
    short int facet_e_count; // value represents count of facet in different layers
};

struct EdgeData
{
    std::vector<bool> edge_e_layers; // vector of edge existence in layers
    size_t id;
    double length;
};

typedef boost::adjacency_list<boost::vecS, boost::hash_setS,
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
typedef boost::graph_traits<BoostFacetGraph>::out_edge_iterator out_edge_iterator;


typedef boost::filtered_graph<BoostFacetGraph, boost::function<bool(BoostFacetGraph::edge_descriptor)>, boost::function<bool(BoostFacetGraph::vertex_descriptor)> > ComponentGraph;
typedef boost::shared_ptr<std::vector<unsigned long>> vertex_component_map;

typedef std::map<vertex_descriptor, size_t> vertex_descriptor_map;