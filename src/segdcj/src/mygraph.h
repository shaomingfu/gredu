#ifndef __MYGRAPH__
#define __MYGRAPH__

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/max_cardinality_matching.hpp>

#include <vector>

using namespace std;
using namespace boost;

// undirected graph
typedef adjacency_list<setS, vecS, undirectedS> ugraph;

typedef graph_traits<ugraph>::vertex_iterator vertex_iterator;
typedef graph_traits<ugraph>::vertex_descriptor vertex_descriptor;
typedef graph_traits<ugraph>::edge_iterator edge_iterator;
typedef graph_traits<ugraph>::out_edge_iterator out_edge_iterator;
typedef graph_traits<ugraph>::edge_descriptor edge_descriptor;

static vertex_descriptor VNULL = graph_traits<ugraph>::null_vertex();

typedef edmonds_augmenting_path_finder<ugraph, int*, int*> path_finder;

// directed graph
typedef adjacency_list<setS, vecS, directedS> dgraph;

typedef graph_traits<dgraph>::vertex_iterator vertex_iterator_d;
typedef graph_traits<dgraph>::vertex_descriptor vertex_descriptor_d;
typedef graph_traits<dgraph>::edge_iterator edge_iterator_d;
typedef graph_traits<dgraph>::out_edge_iterator out_edge_iterator_d;
typedef graph_traits<dgraph>::edge_descriptor edge_descriptor_d;

//typedef property_map<dgraph, vertex_index_t>::const_type const_vertex_index_map_d;

int calc_shortest_cycle(int start, vector<int> & cycle, const dgraph & gr, const vector<bool> & del);
int calc_cycles(vector< vector<int> > & cycles, const dgraph & gr);

#endif
