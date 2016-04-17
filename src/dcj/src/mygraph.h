#ifndef __MYGRAPH__
#define __MYGRAPH__

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/max_cardinality_matching.hpp>

using namespace boost;

typedef property<vertex_index_t, int> vertex_properties;
//typedef property<edge_index_t, int> edge_properties;
typedef adjacency_list<setS, vecS, undirectedS, vertex_properties, no_property> ugraph;

typedef graph_traits<ugraph>::vertex_iterator vertex_iterator;
typedef graph_traits<ugraph>::vertex_descriptor vertex_descriptor;
typedef graph_traits<ugraph>::edge_iterator edge_iterator;
typedef graph_traits<ugraph>::out_edge_iterator out_edge_iterator;
typedef graph_traits<ugraph>::edge_descriptor edge_descriptor;

typedef property_map<ugraph, vertex_index_t>::const_type const_vertex_index_map;
//typedef property_map<ugraph, edge_index_t>::type edge_index_map;

static vertex_descriptor VNULL = graph_traits<ugraph>::null_vertex();

typedef edmonds_augmenting_path_finder<ugraph, vertex_descriptor*, const_vertex_index_map> path_finder;
typedef path_finder::vertex_to_vertex_map_t mate_type;

#endif
