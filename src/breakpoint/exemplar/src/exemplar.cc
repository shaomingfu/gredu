#include "exemplar.h"
#include "psolver.h"
#include "ilp.h"
#include "spliter.h"
#include "decomposer.h"
#include "bbsearch.h"
#include "adjacency.h"

#include <boost/graph/connected_components.hpp>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <algorithm>

using namespace std;

exemplar::exemplar(config * _conf, genome * g1, genome * g2)
	: conf(_conf), gm1(g1), gm2(g2)
{
	alphabet_size = conf->alphabet_size;
	nshadj = 0;
}

exemplar::~exemplar()
{
	alphabet_size = conf->alphabet_size;
}

int exemplar::solve()
{
	set_num_exemplars();

	if(conf->algo == 1) simplify_ilp();
	else if(conf->algo == 2) branch_bound();
	
	return 0;
}

int exemplar::branch_bound()
{
	bbsearch bbs(gm1, gm2);
	bbs.solve();
	return 0;
}

int exemplar::simplify_ilp()
{
	while(true)
	{
		bool b = resolve();
		if(b == false) break;
	}

	relabel(gm1, x2y);
	relabel(gm2, y2x);

	printf("genomes' size = (%d,%d), mapping has %d pairs, %d shared adjacencies\n", gm1->size(), gm2->size(), (int)x2y.size(), gm1->sadist(*gm2, x2y));

	improve();

	printf("genomes' size = (%d,%d), mapping has %d pairs, %d shared adjacencies\n", gm1->size(), gm2->size(), (int)x2y.size(), gm1->sadist(*gm2, x2y));

	return 0;
}

bool exemplar::resolve()
{
	//remove_tandem();
	statistic();

	vector<int> v = nex;

	while(true)
	{
		psolver p(gm1, gm2, conf->max_combination);
		bool b0 = p.extend();
		purify(p.ex1, p.ex2);
		statistic();
		if(!b0) break;
	}

	while(true)
	{
		psolver p(gm1, gm2, conf->max_combination);
		bool b0 = p.extend();
		bool b1 = p.simplify1();
		purify(p.ex1, p.ex2);
		statistic();
		if(!b1 && !b0) break;
	}

	if(conf->heuristic == true) return false;

	while(true)
	{
		psolver p(gm1, gm2, conf->max_combination);
		bool b0 = p.extend();
		bool b1 = p.simplify1();
		bool b2 = p.simplify2();
		purify(p.ex1, p.ex2);
		//bool b3 = p.simplify3();
		//purify_su(gm1, p.su1);
		//purify_su(gm2, p.su2);
		//purify_shared_adjacencies();
		//remove_tandem();
		statistic();
		if(!b0 && !b1 && !b2) break;
	}

	vector<bool> m;
	m.assign(conf->alphabet_size + 1, false);
	for(int i = 1; i <= alphabet_size; i++)
	{
		if(v[i] > nex[i]) m[i] = true;
	}

	vector<PG> p1, p2;
	vector<gene*> d1, d2;
	shrink(gm1, m, p1, d1);
	shrink(gm2, m, p2, d2);

	statistic("SSS");

	ilp * lp = new ilp(gm1, gm2);
	lp->model->getEnv().set(GRB_DoubleParam_TimeLimit, conf->ilp_timelimit);
	lp->model->getEnv().set(GRB_IntParam_MIPFocus, conf->ilp_focus);
	lp->solve();
	printf("\n");

	recover(gm1, p1, d1);
	recover(gm2, p2, d2);

	bool b = purify(lp->ex1, lp->ex2);
	//purify_shared_adjacencies();
	//remove_tandem();
	statistic();
	delete lp;

	return b;
}

int exemplar::set_num_exemplars()
{
	vector<int> v1 = gm1->build_gene_copy();
	vector<int> v2 = gm2->build_gene_copy();

	assert(v1.size() == v2.size());
	nex.assign(v1.size(), 1);

	for(int i = 1; i < v1.size(); i++)
	{
		int n = v1[i] < v2[i] ? v1[i] : v2[i];
		if(conf->num_exemplars >= 0 && conf->num_exemplars < n) n = conf->num_exemplars;
		nex[i] = n;
	}
	return 0;
}

int exemplar::relabel(genome * gm, const MPG & xy)
{
	vector< vector<gene*> > gf;
	gm->build_gene_map(gf);

	for(int f = 1; f < gf.size(); f++)
	{
		for(int i = 0; i < gf[f].size(); i++)
		{
			gene * x = gf[f][i];
			if(fm.find(f) != fm.end())
			{
				if(x->x > 0) x->x = fm[f];
				else x->x = 0 - fm[f];
			}
			if(xy.find(x) == xy.end())
			{
				operation * op = new deletion(x, x, true);
				gm->do_deletion(op);
				delete op;
			}
		}
	}
	return 0;
}

int exemplar::statistic(string s)
{
	int f = 0;
	vector<int> v1 = gm1->build_gene_copy();
	vector<int> v2 = gm2->build_gene_copy();
	for(int i = 1; i <= alphabet_size; i++)
	{
		if(nex[i] == 0) f++;
		//if(nex[i] >= 1) printf("family %5d, nex = %5d, gm1 = %5d, gm2 = %5d\n", i, nex[i], v1[i], v2[i]);
	}
	
	printf("%sgenome sizes = (%6d, %6d), %6d pairs are fixed, %6d out of %6d families are fixed\n\n",
			s.c_str(), gm1->size(), gm2->size(), (int)x2y.size(), f, alphabet_size);

	/*
	gm1->print();
	printf("---\n");
	gm2->print();
	printf("\n");

	printf("================\n");
	gm1->printc(alphabet_size);
	printf("---\n");
	gm2->printc(alphabet_size);
	printf("\n");
	*/

	return 0;
}

bool exemplar::purify(const vector<gene*> & x, const vector<gene*> & y)
{
	assert(x.size() == y.size());

	bool result = false;

	vector< vector<gene*> > gf1;
	vector< vector<gene*> > gf2;
	gm1->build_gene_map(gf1);
	gm2->build_gene_map(gf2);

	for(int i = 1; i <= alphabet_size; i++)
	{
		if(x[i] == NULL) continue;
		assert(y[i] != NULL);

		if(nex[i] == 0) continue;

		gene * gx = x[i];
		gene * gy = y[i];

		//printf("fix %5d, %5d, family = %d, alphabet = %d\n", gx->x, gy->x, nex[i], conf->alphabet_size);

		conf->alphabet_size++;
		nex.push_back(1);
		assert(nex.size() - 1 == conf->alphabet_size);

		fm.insert(PI(conf->alphabet_size, (int)fabs(gx->x)));

		if(gx->x > 0) gx->x = conf->alphabet_size;
		else gx->x = 0 - conf->alphabet_size;

		if(gy->x > 0) gy->x = conf->alphabet_size;
		else gy->x = 0 - conf->alphabet_size;

		assert(x2y.find(gx) == x2y.end());
		assert(y2x.find(gy) == y2x.end());
		x2y.insert(PG(gx, gy));
		y2x.insert(PG(gy, gx));

		result = true;
		nex[i]--;
		
		if(nex[i] >= 1) continue;

		purify_family(gm1, gf1[i], x2y);
		purify_family(gm2, gf2[i], y2x);
	}

	return result;
}

bool exemplar::purify_su(genome * gm, const set<gene*> & su)
{
	bool b = false;
	vector<int> v = gm->build_gene_copy();
	set<gene*>::const_iterator it;
	for(it = su.begin(); it != su.end(); it++)
	{
		gene * g = (*it);
		int f = (int)fabs(g->x);
		assert(v[f] >= nex[f]);
		if(nex[f] == 0) continue;
		//if(v[f] == nex[f]) continue;
		if(v[f] == 1 && nex[f] == 1) continue;
		operation * op = new deletion(g, g, true);
		gm->do_deletion(op);
		delete op;
		b = true;
		if(v[f] == nex[f]) nex[f]--;
		v[f]--;
	}
	return b;
}

int exemplar::purify_family(genome * gm, const vector<gene*> & v, const MPG & xy)
{
	for(int k = 0; k < v.size(); k++)
	{
		gene * g = v[k];
		assert(g->x != 0);
		if(xy.find(g) != xy.end()) continue;
		operation * op = new deletion(g, g, true);
		gm->do_deletion(op);
		delete op;
	}
	return 0;
}

int exemplar::purify_shared_adjacencies()
{
	/*
	printf("----------------------------------------------------------\n\n");
	printf("BEFORE PURIFY SHARED ADJACENCIES\n");
	pbase pb1(gm1, gm2);
	pb1.print_shared_adjacencies();
	*/

	map<gene*, int> m1 = gm1->build_gene_indices();
	map<gene*, int> m2 = gm2->build_gene_indices();
	set<gene*> s1;
	set<gene*> s2;

	MPG::iterator it;
	for(it = x2y.begin(); it != x2y.end(); it++)
	{
		gene * x = it->first;
		gene * y = it->second;
		if(m1.find(x) == m1.end()) continue;
		if(m2.find(y) == m2.end()) continue;
		gene * xx = x->b;
		gene * yy = NULL;
		if(x->x == y->x) yy = y->b;
		else if(x->x + y->x == 0) yy = y->a;
		else assert(false);
		if(x2y.find(xx) == x2y.end()) continue;
		if(x2y[xx] != yy) continue;
		assert(y2x.find(yy) != y2x.end());
		assert(y2x[yy] == xx);
		if(x->x == y->x && xx->x != yy->x) continue;
		if(x->x + y->x == 0 && xx->x + yy->x != 0) continue;
		if(x->x == 0 || xx->x == 0) continue;
		assert(s1.find(xx) == s1.end());
		assert(s2.find(yy) == s2.end());
		s1.insert(xx);
		s2.insert(yy);
	}

	assert(s1.size() == s2.size());
	nshadj += s1.size();

	for(set<gene*>::iterator it = s1.begin(); it != s1.end(); it++)
	{
		gene * g = (*it);
		operation * op = new deletion(g, g, false);
		gm1->do_deletion(op);
		delete op;
	}

	for(set<gene*>::iterator it = s2.begin(); it != s2.end(); it++)
	{
		gene * g = (*it);
		operation * op = new deletion(g, g, false);
		gm2->do_deletion(op);
		delete op;
	}

	/*
	printf("AFTER  PURIFY SHARED ADJACENCIES\n");
	pbase pb2(gm1, gm2);
	pb2.print_shared_adjacencies();
	printf("----------------------------------------------------------\n\n");
	*/

	return 0;
}

int exemplar::remove_tandem()
{
	gm1->remove_tandem();
	gm2->remove_tandem();

	vector<int> v1 = gm1->build_gene_copy();
	vector<int> v2 = gm2->build_gene_copy();

	for(int i = 1; i < v1.size(); i++)
	{
		int n = v1[i] < v2[i] ? v1[i] : v2[i];
		if(nex[i] > n) nex[i] = n;
	}

	return 0;
}

int exemplar::shrink(genome * gm, const vector<bool> & m, vector<PG> & p, vector<gene*> & d)
{
	for(int i = 0; i < gm->chrms.size(); i++)
	{
		chrm * ch = gm->chrms[i];
		assert(ch->type == LINEAR);
		gene * q = ch->p->b;
		gene * s = NULL;
		while(q != NULL)
		{
			int f = (int)fabs(q->x);
			if(m[f] == true) 
			{
				if(s == NULL) s = q;
			}
			else if(s != NULL)
			{
				p.push_back(PG(s, q->a));
				d.push_back(q);
				operation * op = new deletion(s, q->a, false);
				gm->do_deletion(op);
				delete op;
				s = NULL;
			}
			q = q->b;
		}
	}
	return 0;
}

int exemplar::recover(genome * gm, const vector<PG> & p, const vector<gene*> & d)
{
	assert(p.size() == d.size());
	for(int i = 0; i < p.size(); i++)
	{
		gene * x = p[i].first;
		gene * y = p[i].second;
		x->a = NULL;
		y->b = NULL;
		operation * op = new insertion(d[i], x, y);
		gm->do_insertion(op);
		delete op;
	}
	return 0;
}

int exemplar::inc_shadj(gene * x, gene * y)
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

int exemplar::improve()
{
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
