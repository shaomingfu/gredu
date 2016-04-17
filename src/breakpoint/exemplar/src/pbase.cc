#include "pbase.h"
#include "mygraph.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <algorithm>

using namespace std;

pbase::pbase(genome * g1, genome * g2)
	: gm1(g1), gm2(g2)
{
	init();
	build_gene_indices();
	build_adjacencies(gm1, va1);
	build_adjacencies(gm2, va2);
	build_shared_adjacencies();
	build_span_map();
	fix_singletons();
	build_candidates();

	/*
	print_shared_adjacencies();
	gm1->printi(gi);
	printf("---\n");
	gm2->printi(gi);
	printf("\n");
	*/
}

pbase::~pbase()
{
}

int pbase::init()
{
	gm1->build_gene_map(gf1);
	gm2->build_gene_map(gf2);
	assert(gf1.size() == gf2.size());
	ex1.assign(gf1.size(), NULL);
	ex2.assign(gf2.size(), NULL);
	return 0;
}

PG pbase::sort(gene * x, gene * y)
{
	int xx = gi[x];
	int yy = gi[y];
	if(xx < yy) return PG(x, y);
	else return PG(y, x);
}

int pbase::print_shared_adjacencies()
{
	for(int i = 0; i < vsa.size(); i++)
	{
		shadj sa = vsa[i];
		string s = sa.print_string(gi);
		PG p1 = sort(sa.x1, sa.x2);
		PG p2 = sort(sa.y1, sa.y2);
		int d1 = distance(p1);
		int d2 = distance(p2);
		printf("shadj %6d: %s, length = (%2d, %2d)\n", i, s.c_str(), d1, d2);
	}
	return 0;
}

int pbase::build_gene_indices()
{
	gi.clear();
	ig.clear();
	for(int i = 0; i < gm1->chrms.size(); i++)
	{
		gene * q = gm1->chrms.at(i)->p;
		while(q != NULL)
		{
			ig.push_back(q);
			gi.insert(pair<gene*, int>(q, gi.size()));
			q = q->b;
			if(q == gm1->chrms.at(i)->p) break;
		}
	}

	for(int i = 0; i < gm2->chrms.size(); i++)
	{
		gene * q = gm2->chrms.at(i)->p;
		while(q != NULL)
		{
			ig.push_back(q);
			gi.insert(pair<gene*, int>(q, gi.size()));
			q = q->b;
			if(q == gm2->chrms.at(i)->p) break;
		}
	}
	assert(ig.size() == gi.size());
	return 0;
}

int pbase::build_adjacencies(genome * gm, vector<adjacency> & va)
{
	va.clear();
	vector<int> gc = gm->build_gene_copy();
	for(int i = 0; i < gm->chrms.size(); i++)
	{
		chrm * ch = gm->chrms[i];
		assert(ch->type == LINEAR);
		gene * s = ch->p;
		int gs = (int)fabs(s->x);
		while(s->b != NULL)
		{
			gene * t = s->b;
			int gt = (int)fabs(t->x);
			vector<int> v;
			set<int> ss;
			while(true)
			{
				if(s->x == t->x) break;			// TODO, this is a heuristic
				if((int)fabs(s->x) != (int)fabs(t->x) && ss.find(t->x) == ss.end()) 
				{
					va.push_back(adjacency(PG(s, t)));
				}
				if(t->b == NULL) break;
				if(gc[gt] <= 1) break;
				if(ss.find(t->x) == ss.end()) ss.insert(t->x);
				gc[gt]--;
				v.push_back(gt);
				t = t->b;
				gt = (int)fabs(t->x);
			}
			for(int k = 0; k < v.size(); k++) gc[v[k]]++;
			s = s->b;
			gs = (int)fabs(s->x);
		}
	}
	return 0;
}

int pbase::build_shared_adjacencies()
{
	vsa.clear();
	vct.clear();
	vct.resize(gi.size());

	map< int, vector<int> > m;
	for(int i = 0; i < va2.size(); i++)
	{
		adjacency & a = va2[i];
		//if(a.available() == false) continue;
		int x = fabs(a.e1.g->x) + fabs(a.e2.g->x);
		if(m.find(x) == m.end())
		{
			vector<int> v;
			v.push_back(i);
			m.insert(pair<int, vector<int> >(x, v));
		}
		else
		{
			m.find(x)->second.push_back(i);
		}
	}

	for(int i = 0; i < va1.size(); i++)
	{
		adjacency & a = va1[i];
		//if(a.available() == false) continue;
		int x = fabs(a.e1.g->x) + fabs(a.e2.g->x);
		if(m.find(x) == m.end()) continue;
		vector<int> & v = m.find(x)->second;
		for(int j = 0; j < v.size(); j++)
		{
			int jj = v[j];
			adjacency & b = va2[jj];
			shadj sa;
			if(a.e1.g->x == b.e1.g->x && a.e2.g->x == b.e2.g->x)
			{
				sa = shadj(a.e1.g, a.e2.g, b.e1.g, b.e2.g);
			}
			else if(a.e1.g->x + b.e2.g->x == 0 && a.e2.g->x + b.e1.g->x == 0)
			{
				sa = shadj(a.e1.g, a.e2.g, b.e2.g, b.e1.g);
			}
			else continue;

			//if(consistent(sa) == false) continue;

			assert(a.weak_compare(b) >= 2);

			vsa.push_back(sa);

			vct[gi[sa.x1]].insert(vsa.size() - 1);
			vct[gi[sa.x2]].insert(vsa.size() - 1);
			vct[gi[sa.y1]].insert(vsa.size() - 1);
			vct[gi[sa.y2]].insert(vsa.size() - 1);
		}
	}
	
	return 0;
}

int pbase::build_span_map()
{
	vsp.clear();
	vsp.resize(gi.size());
	for(int i = 0; i < vsa.size(); i++)
	{
		shadj & sa = vsa[i];
		PG p = sort(sa.x1, sa.x2);
		gene * g = p.first->b;
		while(g != p.second)
		{
			assert(gi[g] < vsp.size());
			vsp[gi[g]].insert(i);
			g = g->b;
		}

		p = sort(sa.y1, sa.y2);
		g = p.first->b;
		while(g != p.second)
		{
			assert(gi[g] < vsp.size());
			vsp[gi[g]].insert(i);
			g = g->b;
		}
	}
	return 0;
}

int pbase::build_candidates()
{
	for(int i = 0; i < vsa.size(); i++) build_candidate(vsa[i]);
	for(int i = 0; i < cds.size(); i++) cds[i].jp = jumpable(cds[i]);
	return 0;
}

bool pbase::build_candidate(const shadj & sa)
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

bool pbase::jumpable(const candidate & cd)
{
	int x = gi[cd.xv[0]];
	for(set<int>::iterator it = vct[x].begin(); it != vct[x].end(); it++)
	{
		shadj & sa = vsa[*it];
		if(sa.x1 == cd.xv[0] && sa.y1 == cd.yv[0] && gi[sa.x2] < x) return true;
		if(sa.x2 == cd.xv[0] && sa.y2 == cd.yv[0] && gi[sa.x1] < x) return true;
	}

	int n = cd.size() - 1;
	x = gi[cd.xv[n]];
	for(set<int>::iterator it = vct[x].begin(); it != vct[x].end(); it++)
	{
		shadj & sa = vsa[*it];
		if(sa.x1 == cd.xv[n] && sa.y1 == cd.yv[n] && gi[sa.x2] > x) return true;
		if(sa.x2 == cd.xv[n] && sa.y2 == cd.yv[n] && gi[sa.x1] > x) return true;
	}

	return false;
}

int pbase::fix_singletons()
{
	for(int i = 1; i < gf1.size(); i++)
	{
		if(gf1[i].size() != 1) continue;
		if(gf2[i].size() != 1) continue;
		ex1[i] = gf1[i][0];
		ex2[i] = gf2[i][0];
	}
	return 0;
}
