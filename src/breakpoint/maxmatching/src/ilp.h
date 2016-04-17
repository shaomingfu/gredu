#ifndef __ILP_H__
#define __ILP_H__

#include "pbase.h"
#include "genome.h"
#include "gurobi_c++.h"
#include "mygraph.h"
#include <vector>
#include <set>
#include <map>

using namespace std;

class ilp : public pbase
{
public:
	ilp(config * conf, genome * g1, genome * g2, const MPG & xy, const MPG & yx);
	virtual ~ilp();

public:
	// members for ILP
	vector<GRBVar> gvars;					// gene variables
	vector<GRBVar> avars;					// adjacency variables

	GRBModel * model;
	GRBEnv * env;
	int objval;

	set<int> saopt;

public:
	int solve();

	int add_gene_variables();
	int add_shared_adjacency_variables();

	int add_number_constraints();
	int add_inference_constraints();
	int add_inference_constraints(const candidate & c);
	int add_shared_adjacency_edge_constraints();
	int add_bridge_constraint(gene * p, gene * q, int e);
	int add_coexist_constraints();

	int set_objective();

	int collect();
	int update();
	int analyze();
};

#endif
