#ifndef __ILP0_H__
#define __ILP0_H__

#include "pbase0.h"
#include "genome.h"
#include "gurobi_c++.h"
#include "mygraph.h"
#include <vector>
#include <set>
#include <map>

using namespace std;

// implementation of Augibaud's ILP formulation

class ilp0 : public pbase0
{
public:
	ilp0(config * conf, genome * g1, genome * g2, const MPG & xy, const MPG & yx);
	virtual ~ilp0();

public:
	// members for ILP
	vector<GRBVar> gvars;					// gene variables
	vector<GRBVar> pvars;					// pair variables
	vector<GRBVar> avars;					// adjacency variables

	GRBModel * model;
	GRBEnv * env;
	int objval;

public:
	int solve();

	int add_gene_variables();
	int add_pair_variables();
	int add_shared_adjacency_variables();

	int add_number_constraints();
	int add_pair_gene_constraints();
	int add_fixed_pair_constraints();
	int add_pair_coexist_constraints();
	int add_shared_adjacency_pair_constraints();
	int add_bridge_constraint(gene * p, gene * q, int e);
	int add_inference_constraints();

	int set_objective();
	int collect();
};

#endif
