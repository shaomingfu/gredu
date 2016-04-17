#ifndef __DISTOPTIMIZER_H__
#define __DISTOPTIMIZER_H__ 

#include "mygraph.h"
#include "genome.h"
#include "config.h"

#include <vector>
#include <set>

using namespace std;

typedef pair<PI, int> PII;

class maxmatching
{
public:
	maxmatching(config * _conf, genome * _g1, genome * _g2);
	virtual ~maxmatching();

public:
	config * conf;
	
	genome * gm1;
	genome * gm2;

	map<gene*, gene*> x2y;
	map<gene*, gene*> y2x;

	int nshadj;

public:
	int remove_tandem();
	int solve();
	int simplify_ilp();
	int statistic(string s = "");
	int improve();
	int inc_shadj(gene * x, gene * y);
	int print_mapping();
};

#endif
