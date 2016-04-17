#ifndef _BBSEARCH__H__
#define _BBSEARCH__H__

#include "gurobi_c++.h"
#include "candidate.h"
#include "gfamily.h"
#include "bbstate.h"
#include "isolver.h"
#include <vector>
#include <set>
#include <map>

using namespace std;

class bbsearch : public isolver
{
public:
	bbsearch(genome * g1, genome * g2);
	virtual ~bbsearch();

	vector<gfamily> gfs;
	int mind;

	GRBEnv * env;
	vector<candidate> cds;
	vector< vector<int> > cdv;

public:
	int solve();
	int build_families();
	int sort_families();
	int search();
	int init_bbstate(bbstate & bs);
	int fix_candidate(const candidate & cd, bbstate & bs);
	set<int> add_pair(gene * x, gene * y, bbstate & bs);
	set<int> add_pair(const PI & p, bbstate & bs);
	int dfs(const bbstate & bs, int level);
	int optimize(bbstate & bs, const set<int> & s);
	bool optimize(bbstate & bs, int cd);
	int print();

	int build_candidates();
	bool build_candidate(const shadj & sa);
};

#endif
