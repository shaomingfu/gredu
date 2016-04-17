#ifndef __LPSOLVER2_H__
#define __LPSOLVER2_H__

#include "mygraph.h"
#include "config.h"
#include "gurobi_c++.h"

typedef pair<int, int> PI;

class lpsolver2
{
public:
	lpsolver2(config * _conf, const ugraph & _ar, const vector<PI> & _ge,
			const vector<int> & _gf, const vector<int> & _adj,
			const vector< vector<int> > & _vdup, int flag);
	~lpsolver2();

public:
	config * conf;
	GRBModel * model;
	GRBEnv * env;

public:
	const ugraph & ar;
	const vector<PI> & ge;
	const vector<int> & gf;
	const vector<int> & adj;
	const vector< vector<int> > & vdup;

	vector<int> eg;
	map<edge_descriptor, int> ei;
	map<int, edge_descriptor> ie;

	map<int, int> x2y;
	map<int, int> y2x;

public:
	vector<int> vu;			// upper bounds
	vector<GRBVar> gvars;	// variables indicating whether genes are duplicated
	vector<GRBVar> lvars;	// variables representing the labels for each extremity
	vector<GRBVar> xvars;	// variables indicating whther the label reaches its upper bound
	vector<GRBVar> evars;	// variables indicating whether each edge is chosen
	vector<GRBVar> ovars;	// variables indicating whether each operation is chosen

public:
	double solve();

	// collect returning operations
	int statistic_solution();
	int collect_pairs();
	int collect_fixed_pairs();
	int print_solution();
	int get_num_cycles();

private:
	// build data structures
	int build_gene_extremities();
	int build_edge_map();
	int build_upper_bounds();

	// add variables
	int add_gene_variables();
	int add_label_variables();
	int add_bound_variables();
	int add_edge_variables();
	int add_operation_variables();

	// add constraints
	// TODO, two edges corresponding to the same gene should have the same label
	int add_gene_label_constraints();
	int add_adjacency_constraints();
	int add_edge_label_constraints();
	int add_bound_constraints();
	int add_family_constraints();
	int add_gene_edge_constraints();
	int add_gene_operation_constraints();
	int add_operation_constraints();
	int add_consistent_constraints();
	
	// set objective function
	int set_objective();

};

#endif
