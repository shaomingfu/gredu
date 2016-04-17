#ifndef __PSOLVER_H__
#define __PSOLVER_H__

#include "pbase.h"
#include "gurobi_c++.h"
#include <vector>
#include <set>
#include <map>

using namespace std;

// only for gene size = 1

class psolver : public pbase
{
public:
	psolver(genome * g1, genome * g2, int _max_combination);
	virtual ~psolver();

public:
	set<gene*> su1;
	set<gene*> su2;
	vector<bool> cdj;

	GRBEnv * env;
	int max_combination;

public:
	bool extend();
	bool simplify1();
	bool simplify2();
	bool simplify3();

	bool extend(const shadj & sa);
	bool simplify1(int cdi);
	bool simplify2(int f);
	bool simplify3(const vector<gene*> & gf, set<gene*> & su);
	bool testify(const vector<gene*> & xv, const vector<gene*> & yv);

	bool innocent(const shadj & sa, const vector<gene*> & xv, const vector<gene*> & yv);
	set<int> build_span_intersection(const vector<gene*> & v);
	set<int> build_span_union(const vector<gene*> & v);
	set<int> build_contact_union(const vector<gene*> & v);
	set<int> set_union(const set<int> & x, const set<int> & y);
	set<int> set_intersection(const set<int> & x, const set<int> & y);

	int statistic();
};

#endif
