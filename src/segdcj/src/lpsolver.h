#ifndef __LPSOLVER_H__
#define __LPSOLVER_H__

#include "psolver.h"
#include "mygraph.h"
#include "genome.h"
#include "config.h"
#include "gurobi_c++.h"

class lpsolver : public psolver
{
public:
	lpsolver(config * _conf, genome * _gm1, genome * _gm2);
	~lpsolver();

public:
	GRBModel * model;
	GRBEnv * env;

public:
	// upper bounds
	vector<int> vu1;
	vector<int> vu2;

	// variables indicating whether genes are duplicated
	vector<GRBVar> gvars1;
	vector<GRBVar> gvars2;

	// variables representing the labels for each extremity
	vector<GRBVar> lvars1;
	vector<GRBVar> lvars2;

	// variables indicating whther the label reaches its upper bound
	vector<GRBVar> xvars1;
	vector<GRBVar> xvars2;

	// variables indicating whether each edge is chosen
	vector<GRBVar> evars;

	// variables indicating whether each operation is chosen
	vector<GRBVar> ovars1;
	vector<GRBVar> ovars2;

	map<gene*, gene*> x2y;
	map<gene*, gene*> y2x;

public:
	int solve();

private:
	int build_upper_bounds();

	// add variables
	int add_gene_variables(vector<GRBVar> & gvars, int n);
	int add_label_variables(vector<GRBVar> & lvars, int n, const vector<int> & vu);
	int add_bound_variables(vector<GRBVar> & xvars, int n);
	int add_edge_variables();
	int add_operation_variables(vector<GRBVar> & ovars, int n);

	// add constraints
	int add_gene_label_constraints(const vector<GRBVar> & gvars, const vector<GRBVar> & lvars, const vector<int> & vu);
	int add_adjacency_constraints(map<gene*, int> & gi, const vector<GRBVar> & gvars, const vector<GRBVar> & lvars, const vector<int> & vu);
	int add_edge_label_constraints();
	int add_bound_constraints(const vector<GRBVar> & lvars, const vector<GRBVar> & xvars, const vector<int> & vu);
	int add_family_constraints();
	int add_gene_edge_constraints();
	int add_gene_operation_constraints(map<gene*, int> & gi, const vector<PG> & dup, const vector<GRBVar> & ovars, const vector<GRBVar> & gvars);
	int add_operation_constraints(map<gene*, int> & gi, const vector<PG> & dup, const vector<GRBVar> & ovars, const vector<GRBVar> & gvars, const vector<GRBVar> & lvars, const vector<int> & vu);
	int add_null_gene_constraints(const vector<gene*> & v, map<gene*, int> & gi, const vector<GRBVar> & lvars);
	int add_inference_constraints();
	
	// set objective function
	int set_objective1();
	int set_objective2();

	// collect returning operations
	int statistic_solution();
	int write_duplicated_segments();
	int collect_pairs();
	int collect_fixed_pairs();
	int remove_duplicated_genes();
	int print();
	int print_solution();

};

#endif
