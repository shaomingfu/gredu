#ifndef __ILP1_H__
#define __ILP1_H__

#include "spliter.h"
#include "genome.h"
#include "gurobi_c++.h"
#include "mygraph.h"
#include <vector>
#include <set>
#include <map>

using namespace std;

class ilp : public spliter
{
public:
	ilp(genome * g1, genome * g2);
	virtual ~ilp();

public:
	// members for ILP
	vector<GRBVar> gvars;					// gene variables
	vector<GRBVar> avars;					// adjacency variables

	GRBModel * model;
	GRBEnv * env;
	int objval;

public:
	int solve();

	int add_gene_variables();
	int add_shared_adjacency_variables();

	int add_number_constraints();
	int add_inference_constraints();
	int add_inference_constraints(const candidate & c);
	int add_shared_adjacency_edge_constraints();
	int add_bridge_constraint(gene * p, gene * q, int e);

	int set_objective();
	int collect();
};

#endif
