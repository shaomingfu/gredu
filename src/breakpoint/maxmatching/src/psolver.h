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
	psolver(config * conf, genome * g1, genome * g2, const MPG & xy, const MPG & yx);
	virtual ~psolver();

public:
	set<gene*> su1;
	set<gene*> su2;
	vector<bool> cdj;

	GRBEnv * env;

public:
	bool extend();
	bool extend(const shadj & sa);

	bool simplify1();
	bool simplify1(int cdi);

	bool simplify2();
	bool simplify2(int f);

	bool remove_redundant_genes();

	bool testify(const vector<gene*> & xv, const vector<gene*> & yv);

	int check(const set<int> & saopt);
};

#endif
