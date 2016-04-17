#ifndef __PRESOLVER_H__
#define __PRESOLVER_H__

#include "gurobi_c++.h"

#include "mygraph.h"
#include <vector>
#include <set>
#include <map>

using namespace std;

typedef pair<int, int> PI;
typedef pair<PI, int> PII;

class presolver
{
public:
	presolver(ugraph & _gr, const vector<PII> & _gene_list, const vector<int> & _partners);
	virtual ~presolver();

public:
	ugraph & gr;
	ugraph ig;			// identical graph
	ugraph cg;

	const vector<int> & partners;
	vector<int> nb_indices;
	vector<PI> neighbors;

	const vector<PII> & gene_list;
	vector<int> gn_indices;
	vector<int> complements;

public:
	bool simplify(int level);

	int build_neighbor_indices();
	int build_gene_indices();
	int build_gene_pairs();
	int build_complements();

	int build_identical_graph();
	int analyze_components(const ugraph & g);

	bool simplify_graph_by_fixing();
	bool simplify_graph_by_length2();
	bool simplify_graph_by_half_fixed_length2();
	bool simplify_graph_by_supported_length2(int level);
	bool simplify_graph_by_checked_length2(int level);
	PI calc_identity_degree(int x, int y);
	PI calc_extended_identity_degree(int x, int y);

	int calc_max_identity_degree(int x, int y, const set<int> & xv, const set<int> & yv);
	PI calc_max_identity_degree(int x, int y);

	bool simplify_graph_by_length4();

	bool simplify_graph_by_length6();
	bool fix_length4_cycle(int x);
	bool fix_length6_cycle(int x);

	bool simplify_graph_by_pure_length4();
	bool fix_pure_length4_cycle(int x);
	bool unique_length4_cycle(int x);

	bool simplify_graph_by_components();
	bool fix_component(const vector<int> & v);
	bool remove_extra_edge(int x, const vector<int> & vy);
	int remove_extremity_pair(int x, int y);

	bool is_fixed(int x, int y);

	int fix_adjacency_pair(int x, int y);
	int fix_extremity_pair(int x, int y);
	int fix_extremity(int x);

	int check_identity(int x, int y);				// -1: not identical, 0: identical, 1: double identical

	int print();

	int build_connected_graph();
	bool merge_connected_graph();
	set<int> get_linked_adjacencies(int x);

	int build_subgraph(int x, int y, ugraph & g1, ugraph & g2, vector<PII> & gl, vector<int> & pt);
};

#endif
