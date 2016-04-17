#ifndef __ILP6_H__
#define __ILP6_H__

#include "gurobi_c++.h"
#include "ilp_base.h"

#include "mygraph.h"
#include <vector>
#include <set>
#include <map>

using namespace std;

// single-edge variable, O(n) implementation
class ilp6 : public ilp_base
{
public:
	ilp6(const ugraph & _gr, const vector<PII> & _gene_list, const vector<int> & _partners);
	virtual ~ilp6();

public:
	vector<GRBVar> evars;					// vtype variables
	vector<GRBVar> lvars;					// label variables
	vector<GRBVar> rvars;					// representative variables

	double max_label_index;
	vector<double> ub;

public:
	int solve();

	int set_upper_bounds();
	int add_edge_variables();
	int add_label_variables();
	int add_representative_variables();

	int add_degree_constraints();
	int add_connectivity_constraints();
	int add_label_representative_constraints();
	int add_extra_constraints();

	int set_objective();
	int collect_pairs();
};

#endif
