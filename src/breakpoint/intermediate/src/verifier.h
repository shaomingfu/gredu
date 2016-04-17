#ifndef __VERIFIER_H__
#define __VERIFIER_H__

#include "pbase.h"
#include "gurobi_c++.h"
#include <vector>
#include <set>
#include <map>

using namespace std;

typedef pair<int, int> PI;

class verifier
{
public:
	verifier(const vector<gene*> & x, const vector<gene*> & v, bool dd, pbase * pp, GRBEnv * e, int max_combination);
	virtual ~verifier();

public:
	const vector<gene*> & xv;
	const vector<gene*> & yv;
	bool d;
	pbase * p;
	int max_combination;

	GRBEnv * env;

public:
	bool feasible;			// whether the model is two big
	vector<int> v0;			// given shared adjacencies
	set<int> ssa;			// set of all related shared adjacencies
	map<int, int> sam1;		// innocent shared adjacencies
	map<int, int> sam2;		// affected shared adjacencies
	vector<int> sav1;
	vector<int> sav2;

	int ans1;
	int ans2;
	vector<int> opt1;
	vector<int> opt2;

	map<gene*, int> mi;		// affected genes

public:
	bool verify();
	int build_related_shared_adjacencies();
	int build_affected_shared_adjacencies();
	int build_affected_genes();
	int build_innocent_shared_adjacencies();

	int optimize_innocent_shared_adjacencies();
	int optimize_affected_shared_adjacencies();

	int print_optimal_solutions();
	
	int draw(const string & file);
};

 #endif
