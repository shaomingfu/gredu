#ifndef __DISTOPTIMIZER_H__
#define __DISTOPTIMIZER_H__ 

#include "mygraph.h"
#include "genome.h"
#include "config.h"

#include <vector>
#include <set>

using namespace std;

typedef pair<PI, int> PII;

class exemplar
{
public:
	exemplar(config * _conf, genome * _g1, genome * _g2);
	virtual ~exemplar();

public:
	config * conf;
	
	genome * gm1;
	genome * gm2;

	int alphabet_size;
	vector<int> nex;
	
	int nshadj;

	map<gene*, gene*> x2y;
	map<gene*, gene*> y2x;

	map<int, int> fm;

public:
	int solve();
	int simplify_ilp();
	int branch_bound();
	bool resolve();
	int set_num_exemplars();
	int relabel(genome * gm, const MPG & xy);
	bool purify(const vector<gene*> & x, const vector<gene*> & y);
	bool purify_su(genome * gm, const set<gene*> & s);
	int purify_family(genome * gm, const vector<gene*> & v, const MPG & xy);
	int purify_shared_adjacencies();
	int remove_tandem();
	int statistic(string s = "");
	int shrink(genome * gm, const vector<bool> & m, vector<PG> & p, vector<gene*> & d);
	int recover(genome * gm, const vector<PG> & p, const vector<gene*> & d);
	int inc_shadj(gene * x, gene * y);
	int improve();
};

#endif
