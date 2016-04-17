#include "maxmatching.h"
#include "psolver.h"
#include "ilp.h"
#include "ilp0.h"
#include "adjacency.h"

#include <boost/graph/connected_components.hpp>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <algorithm>

using namespace std;

maxmatching::maxmatching(config * _conf, genome * g1, genome * g2)
	: conf(_conf), gm1(g1), gm2(g2)
{
}

maxmatching::~maxmatching()
{
}

int maxmatching::solve()
{
	simplify_ilp();
	return 0;
}

int maxmatching::simplify_ilp()
{
	//remove_tandem();
	//statistic();

	while(conf->use_psolver)
	{
		psolver p(conf, gm1, gm2, x2y, y2x);
		bool b0 = p.extend();
		bool b1 = false;		//p.simplify1();
		bool b2 = false;		//p.simplify2();

		x2y = p.x2y;
		y2x = p.y2x;
		statistic();
		if(!b0 && !b1 && !b2) break;
	}

	/*
	pbase pb(gm1, gm2, x2y, y2x);
	gm1->printv(pb.gi);
	gm2->printv(pb.gi);
	pb.print_shared_adjacencies();
	pb.print_families();
	*/

	if(conf->heuristic) return 0;

	if(conf->algo == 1)
	{
		ilp * lp = new ilp(conf, gm1, gm2, x2y, y2x);
		lp->solve();
		x2y = lp->x2y;
		y2x = lp->y2x;
		nshadj = lp->objval;
		delete lp;
		printf("\n");
	}
	else if(conf->algo == 0)
	{
		ilp0 * lp = new ilp0(conf, gm1, gm2, x2y, y2x);
		lp->solve();
		x2y = lp->x2y;
		y2x = lp->y2x;
		nshadj = lp->objval;
		delete lp;
		printf("\n");
	}

	printf("genomes size = (%d,%d), %d pairs, (%d,%d) shared adjacencies\n", 
			gm1->size(), gm2->size(), (int)x2y.size(), nshadj, gm1->sadist(*gm2, x2y));

	return 0;
}


int maxmatching::statistic(string s)
{
	int f = 0;
	vector<int> v1 = gm1->build_gene_copy();
	vector<int> v2 = gm2->build_gene_copy();
	
	printf("%sgenome sizes = (%6d, %6d), %6d pairs are fixed\n\n",
			s.c_str(), gm1->size(), gm2->size(), (int)x2y.size());

	return 0;
}

int maxmatching::remove_tandem()
{
	gm1->remove_tandem();
	gm2->remove_tandem();

	return 0;
}

int maxmatching::inc_shadj(gene * x, gene * y)
{
	assert(x->x != 0 && y->x != 0);
	gene * x1 = x->a;
	gene * x2 = x->b;
	gene * y1 = y->a;
	gene * y2 = y->b;
	if(x1->x == 0 || x2->x == 0) return 0;
	if(y1->x == 0 || y2->x == 0) return 0;

	int ans = 0;
	if(x->x == y->x)
	{
		if(x1->x == y1->x) ans--;
		if(x2->x == y2->x) ans--;
	}
	else if(x->x + y->x == 0)
	{
		if(x1->x + y2->x == 0) ans--;
		if(x2->x + y1->x == 0) ans--;
	}
	else assert(false);

	if(adjacency(PG(x1, x2)).weak_compare(adjacency(PG(y1, y2))) >= 2) ans++;

	assert(x2y.find(x1) != x2y.end());
	assert(x2y.find(x2) != x2y.end());

	gene * xx1 = x2y[x1];
	gene * xx2 = x2y[x2];

	if(xx1->b == xx2 && adjacency(PG(x1, x2)).weak_compare(adjacency(PG(xx1, xx2))) >= 2) ans++;
	if(xx2->b == xx1 && adjacency(PG(x1, x2)).weak_compare(adjacency(PG(xx2, xx1))) >= 2) ans++;

	assert(y2x.find(y1) != y2x.end());
	assert(y2x.find(y2) != y2x.end());

	gene * yy1 = y2x[y1];
	gene * yy2 = y2x[y2];

	if(yy1->b == yy2 && adjacency(PG(y1, y2)).weak_compare(adjacency(PG(yy1, yy2))) >= 2) ans++;
	if(yy2->b == yy1 && adjacency(PG(y1, y2)).weak_compare(adjacency(PG(yy2, yy1))) >= 2) ans++;

	return ans;
}

int maxmatching::improve()
{
	// local improving through trying removing every pair of fixed genes
	while(true)
	{
		bool b = false;
		map<gene*, gene*>::iterator it;
		for(it = x2y.begin(); it != x2y.end(); it++)
		{
			gene * x = it->first;
			gene * y = it->second;
			assert(y2x.find(y) != y2x.end());
			assert(y2x[y] == x);

			int inc = inc_shadj(x, y);

			if(inc <= 0) continue;
			b = true;
			operation * opx = new deletion(x, x, true);
			operation * opy = new deletion(y, y, true);
			gm1->do_deletion(opx);
			gm2->do_deletion(opy);
			x2y.erase(x);
			y2x.erase(y);
			delete opx;
			delete opy;
			break;
		}
		if(b == false) break;
	}
	return 0;
}

int maxmatching::print_mapping()
{
	printf("pairs in x2y and y2x:\n");
	MPG::iterator it;
	for(it = x2y.begin(); it != x2y.end(); it++)
	{
		gene * x = it->first;
		gene * y = it->second;
		assert(y2x.find(y) != y2x.end());
		assert(y2x[y] == x);
		printf("genes = (%4d,%4d), names = (%5s,%5s)\n", x->x, y->x, x->s.c_str(), y->s.c_str());
	}
	return 0;
}

