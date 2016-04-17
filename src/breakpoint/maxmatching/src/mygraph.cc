#include "mygraph.h"
#include <cstdio>
#include <iostream>
#include <boost/graph/max_cardinality_matching.hpp>

using namespace std;

int test_graph()
{
	int n = 5;
	ugraph gr(5);

	add_edge(0, 2, gr);
	add_edge(0, 4, gr);
	add_edge(1, 3, gr);
	add_edge(1, 4, gr);

	vertex_iterator vi1, vi2;
	for(tie(vi1, vi2) = vertices(gr); vi1 != vi2; vi1++)
	{
		cout<<*vi1<<endl;
	}

	edge_iterator ei1, ei2;
	for(tie(ei1, ei2) = edges(gr); ei1 != ei2; ei1++)
	{
		cout<<*ei1<<endl;
	}
	
	vector<int> v;
	v.assign(num_vertices(gr), -1);
	v[1] = 4;
	v[4] = 1;

	typedef edmonds_augmenting_path_finder<ugraph, int*, const_vertex_index_map> path_finder;
	path_finder finder(gr, &v[0], get(vertex_index, gr));
	bool b = finder.augment_matching();
	
	printf("b = %d\n", b ? 1 : 0);

	vector<vertex_descriptor> vv(num_vertices(gr));
	property_map<ugraph, vertex_index_t>::type xxx = get(vertex_index, gr);

	path_finder::vertex_to_vertex_map_t mt(vv.begin(), xxx);

	finder.get_current_matching(mt);

	for(int i = 0; i < num_vertices(gr); i++)
	{
		cout<<i<<" -> "<<mt[i]<<endl;
	}

	return 0;
}
