#ifndef __ISOLVER_H__
#define __ISOLVER_H__

#include "pbase.h"
#include "gurobi_c++.h"
#include <vector>
#include <set>
#include <map>

using namespace std;

// only for gene size = 1

class isolver : public pbase
{
public:
	isolver(genome * g1, genome * g2);
	virtual ~isolver();

public:
	GRBEnv * env;

public:
	bool testify(const vector<gene*> & xv, const vector<gene*> & yv, vector<int> & sol, const vector<bool> & sab);
	bool innocent(const shadj & sa, const vector<gene*> & xv, const vector<gene*> & yv);
	set<int> build_span_intersection(const vector<gene*> & v);
	set<int> build_span_union(const vector<gene*> & v);
	set<int> build_contact_union(const vector<gene*> & v);
	set<int> set_union(const set<int> & x, const set<int> & y);
	set<int> set_intersection(const set<int> & x, const set<int> & y);
};

#endif
