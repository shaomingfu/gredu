#include "bbsearch.h"
#include <algorithm>
#include <cmath>

bbsearch::bbsearch(genome * gm1, genome * gm2)
	: isolver(gm1, gm2), mind(9999999)
{
	env = new GRBEnv();
}

bbsearch::~bbsearch()
{
	delete env;
}

int bbsearch::solve()
{
	build_families();
	sort_families();
	build_candidates();
	search();
	return 0;
}

int bbsearch::build_families()
{
	assert(gf1.size() == gf2.size());
	// TODO ignore gene family 0
	for(int i = 1; i < gf1.size(); i++)
	{
		gfamily gf(i);
		for(int j = 0; j < gf1[i].size(); j++)
		{
			int x = gi[gf1[i][j]];
			for(int k = 0; k < gf2[i].size(); k++)
			{
				int y = gi[gf2[i][k]];
				gf.vp.push_back(PI(x, y));
			}
		}
		gfs.push_back(gf);
	}
	return 0;
}

int bbsearch::sort_families()
{
	std::sort(gfs.begin(), gfs.end());
	return 0;
}

int bbsearch::search()
{
	bbstate bs(ig);
	init_bbstate(bs);
	bs.print();
	dfs(bs, 0);
	return 0;
}

int bbsearch::dfs(const bbstate & bs, int level)
{
	if(bs.sx.size() >= gfs.size())
	{
		if(bs.distance < mind) mind = bs.distance;
		printf("*** ");
		bs.print();
		return 0;
	}

	int fi = -1;
	for(int i = 0; i < gfs.size(); i++)
	{
		int f = gfs[i].f;
		if(bs.fx[f] != -1) continue;
		assert(bs.fy[f] == -1);
		fi = i;
		break;
	}

	assert(fi != -1);

	gfamily & gf = gfs[fi];

	for(int k = 0; k < gf.vp.size(); k++)
	{
		gf.vp[k].s = bs.inc_distance(gf.vp[k].p);
	}
	gf.sort();

	for(int k = 0; k < gf.vp.size(); k++)
	{
		if(bs.distance + gf.vp[k].s >= mind) continue;
		bbstate bs2(bs);
		set<int> s = add_pair(gf.vp[k].p, bs2);
		optimize(bs2, s);

		printf("level = %5d, ", level);
		bs2.print();
		printf("\n");

		dfs(bs2, level + 1);
	}

	return 0;
}

int bbsearch::optimize(bbstate & bs, const set<int> & s)
{
	set<int> ss;
	int num = 0;
	int opt = 0;
	for(set<int>::iterator it = s.begin(); it != s.end(); it++)
	{
		int sa = *it;
		set<int> & scd = bs.cdm[sa];
		for(set<int>::iterator ix = scd.begin(); ix != scd.end(); ix++)
		{
			int cd = *ix;
			if(ss.find(cd) != ss.end()) continue;
			ss.insert(cd);
			if(bs.cdb[cd] == false) continue;
			bool b = optimize(bs, cd);
			if(b == true) opt++;
			num++;
		}
	}
	printf("%5d / %5d tests succeed\n", opt, num);
	return 0;
}

bool bbsearch::optimize(bbstate & bs, int cd)
{
	vector<int> sol;
	bool b = testify(cds[cd].xv, cds[cd].yv, sol, bs.sab);
	if(b == true)
	{
		fix_candidate(cds[cd], bs);
	}
	else
	{
		for(int k = 0; k < sol.size(); k++)
		{
			int sa = sol[k];
			bs.affect_candidate(sa, cd);
		}
	}
	return b;
}

int bbsearch::init_bbstate(bbstate & bs)
{
	bs.fx.assign(gf1.size(), -1);
	bs.fy.assign(gf1.size(), -1);

	bs.cdb.resize(cds.size(), true);
	bs.sab.assign(vsa.size(), true);

	for(int i = 0; i < cds.size(); i++)
	{
		if(bs.cdb[i] == false) continue;
		optimize(bs, i);
	}
	return 0;
}

int bbsearch::fix_candidate(const candidate & cd, bbstate & bs)
{
	for(int i = 0; i < cd.xv.size(); i++)
	{
		add_pair(cd.xv[i], cd.yv[i], bs);
	}
	return 0;
}

set<int> bbsearch::add_pair(const PI & p, bbstate & bs)
{
	return add_pair(ig[p.first], ig[p.second], bs);
}

set<int> bbsearch::add_pair(gene * gx, gene * gy, bbstate & bs)
{
	int x = gi[gx];
	int y = gi[gy];

	bs.add_pair(PI(x, y));

	int f = (int)fabs(gx->x);

	// mark unusable candidates
	for(int i = 0; i < cdv[f].size(); i++)
	{
		int c = cdv[f][i];
		bs.cdb[c] = false;
	}
	
	set<int> s;
	//mark unusable shared adjacencies
	for(set<int>::iterator it = vsp[x].begin(); it != vsp[x].end(); it++)
	{
		if(bs.sab[*it] == false) continue;
		s.insert(*it);
		bs.sab[*it] = false;
	}

	for(set<int>::iterator it = vct[x].begin(); it != vct[x].end(); it++)
	{
		shadj & sa = vsa[*it];
		if(bs.sab[*it] == false) continue;
		if(sa.x1 == gx && sa.y1 != gy) 
		{
			s.insert(*it);
			bs.sab[*it] = false;
		}
		else if(sa.x2 == gx && sa.y2 != gy)
		{
			s.insert(*it);
			bs.sab[*it] = false;
		}
	}

	for(set<int>::iterator it = vct[y].begin(); it != vct[y].end(); it++)
	{
		shadj & sa = vsa[*it];
		if(bs.sab[*it] == false) continue;
		if(sa.x1 != gx && sa.y1 == gy) 
		{
			s.insert(*it);
			bs.sab[*it] = false;
		}
		else if(sa.x2 != gx && sa.y2 == gy)
		{
			s.insert(*it);
			bs.sab[*it] = false;
		}
	}
	
	return s;
}

int bbsearch::build_candidates()
{
	for(int i = 0; i < vsa.size(); i++) build_candidate(vsa[i]);

	cdv.resize(gf1.size());
	for(int i = 0; i < cds.size(); i++)
	{
		for(int k = 0; k < cds[i].xv.size(); k++)
		{
			int f = (int)fabs(cds[i].xv[k]->x);
			cdv[f].push_back(i);
		}
	}

	return 0;
}

bool bbsearch::build_candidate(const shadj & sa)
{
	if(sa.adjacent() == false) return false;

	int f1 = (int)fabs(sa.x1->x);
	int f2 = (int)fabs(sa.x2->x);

	if(f1 == f2) return false;

	if(ex1[f1] != NULL && ex1[f2] != NULL) return false;
	assert(ex2[f1] == NULL || ex2[f2] == NULL);

	if(ex1[f1] != NULL && ex1[f1] != sa.x1) return false;
	if(ex1[f2] != NULL && ex1[f2] != sa.x2) return false;
	if(ex2[f1] != NULL && ex2[f1] != sa.y1) return false;
	if(ex2[f2] != NULL && ex2[f2] != sa.y2) return false;

	if(ex1[f1] == sa.x1 && ex2[f1] == sa.y1) return false;
	if(ex1[f2] == sa.x2 && ex2[f2] == sa.y2) return false;

	bool b = sa.direction();
	gene * x = sa.x1;
	gene * y = sa.y1;

	if(b == true && x->a != NULL && y->a != NULL && x->a->x == y->a->x && ex1[(int)fabs(x->a->x)] == NULL) return false;
	if(b == false && x->a != NULL && y->b != NULL && x->a->x + y->b->x == 0 && ex1[(int)fabs(x->a->x)] == NULL) return false;

	vector<gene*> xv;
	vector<gene*> yv;

	vector<bool> m;
	m.assign(gf1.size(), false);
	while(true)
	{
		if(m[(int)fabs(x->x)] == true) break;
		m[(int)fabs(x->x)] = true;
		xv.push_back(x);
		yv.push_back(y);
		x = x->b;
		if(b) y = y->b;
		else y = y->a;

		if(x == NULL || y == NULL) break;
		if(b == true && x->x != y->x) break;
		if(b == false && x->x + y->x != 0) break;
		int f = (int)fabs(x->x);
		if(ex1[f] != NULL) break;
	}

	candidate cd(b, xv, yv);
	cds.push_back(cd);

	return true;
}

int bbsearch::print()
{
	for(int i = 0; i < gfs.size(); i++)
	{
		printf("gene family %5d has %5d pairs\n", gfs[i].f, gfs[i].size());
	}
	printf("\n");

	for(int i = 0; i < cds.size(); i++)
	{
		printf("%5d:", i + 1);
		cds[i].print();	
	}

	return 0;
}

