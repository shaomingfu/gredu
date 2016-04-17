#ifndef __ILP_BASE_H__
#define __ILP_BASE_H__

#include "gurobi_c++.h"

#include "mygraph.h"
#include <vector>
#include <set>
#include <map>

using namespace std;

typedef pair<int, int> PI;
typedef pair<PI, int> PII;

class ilp_base
{
public:
	ilp_base(const ugraph & _gr, const vector<PII> & _gene_list, const vector<int> & _partners);
	virtual ~ilp_base();

public:
	const ugraph & gr;
	const vector<PII> & gene_list;
	const vector<int> & partners;

	vector<int> nb_indices;
	vector<PI> neighbors;

	vector<PI> gene_pairs;
	vector<int> gn_indices;
	vector<int> complements;

	map<edge_descriptor, int> eim;

	GRBModel * model;
	GRBEnv * env;

	// parameters for gurobi solver
	int numericfocus;
	int scaleflag;
	int nodes;
	double gap;
	int cuts;
	int prepasses;
	int gomorypasses;
	int mipfocus;
	double heuristics;

	// result
	set<PI> pairs;

public:
	virtual int solve() = 0;
	int verify_bipartite_graph();
	int build_neighbor_indices();

	int build_complements();
	int build_gene_indices();
	int build_gene_pairs();
	int build_edge_map();

	int set_timelimit(double time);
	int set_parameters();

	int check_connectivity();

	virtual int print();
};

class pedge
{
public:
	pedge(int _a, int _b, int _c, int _d) { a = _a; b = _b; c = _c; d = _d; }
	pedge(const pedge & p) { a = p.a; b = p.b; c = p.c; d = p.d; }
	int print() const { printf("(%5d,%5d,%5d,%5d)", a, b, c, d); return 0; }

public:
	int a, b, c, d;
};

#endif
