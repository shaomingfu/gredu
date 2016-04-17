#include "pbase.h"
#include "mygraph.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <string>

using namespace std;

pbase::pbase(config * c, genome * g1, genome * g2, const MPG & _x2y, const MPG & _y2x)
	: conf(c), gm1(g1), gm2(g2), x2y(_x2y), y2x(_y2x)
{
	init();
}

pbase::~pbase()
{
}

int pbase::init()
{
	gm1->build_gene_map(gf1);
	gm2->build_gene_map(gf2);
	assert(gf1.size() == gf2.size());

	init_fixed_families();
	build_gene_indices();
	build_adjacencies(gm1, va1);
	build_adjacencies(gm2, va2);
	build_shared_adjacencies();
	build_contact_map();

	if(conf->remove_redundant == true)
	{
		remove_redundant_shared_adjacencies();
		build_contact_map();
	}

	build_span_map();
	fix_singletons();
	build_candidates();

	return 0;
}

PG pbase::sort(gene * x, gene * y)
{
	int xx = gi[x];
	int yy = gi[y];
	if(xx < yy) return PG(x, y);
	else return PG(y, x);
}

int pbase::init_fixed_families()
{
	vf.resize(gf1.size(), 0);
	for(MPG::iterator it = x2y.begin(); it != x2y.end(); it++)
	{
		gene * x = it->first;
		int f = (int)fabs(x->x);
		vf[f]++;
	}
	return 0;
}

int pbase::print_families()
{
	for(int i = 0; i < vf.size(); i++)
	{
		int n1 = 0;
		int n2 = 0;
		for(int k = 0; k < gf1[i].size(); k++)
		{
			gene * g = gf1[i][k];
			int x = vct[gi[g]].size();
			if(x >= 1) n1++;
		}
		for(int k = 0; k < gf2[i].size(); k++)
		{
			gene * g = gf2[i][k];
			int x = vct[gi[g]].size();
			if(x >= 1) n2++;
		}

		printf("family %6d ( %3d %3d ) pairs, %2d pairs are fixed, ( %2d %2d ) are in shadj\n",
				i, (int)gf1[i].size(), (int)gf2[i].size(), vf[i], n1, n2);
	}
	return 0;
}

int pbase::print_shared_adjacencies(const set<int> & s)
{
	for(set<int>::const_iterator it = s.begin(); it != s.end(); it++)
	{
		int i = *it;
		shadj sa = vsa[i];
		string s = sa.print_string(gi);
		PG p1 = sort(sa.x1, sa.x2);
		PG p2 = sort(sa.y1, sa.y2);
		int d1 = distance(p1);
		int d2 = distance(p2);

		printf("shadj %4d: %s, length = (%2d, %2d)\n", i, s.c_str(), d1, d2);

		int cx1 = vct[gi[sa.x1]].size();
		int cx2 = vct[gi[sa.x2]].size();
		int cy1 = vct[gi[sa.y1]].size();
		int cy2 = vct[gi[sa.y2]].size();

		int cxy1 = set_intersection(vct[gi[sa.x1]], vct[gi[sa.y1]]).size();
		int cxy2 = set_intersection(vct[gi[sa.x2]], vct[gi[sa.y2]]).size();

		//printf("shadj %4d: %s, length = (%2d, %2d), contact = (%2d,%2d), (%2d,%2d), contact-pairs = (%2d, %2d)\n", 
		//		i, s.c_str(), d1, d2, cx1, cx2, cy1, cy2, cxy1, cxy2);
	}
	return 0;
}

int pbase::print_shared_adjacencies()
{
	set<int> s;
	for(int i = 0; i < vsa.size(); i++) s.insert(i);
	print_shared_adjacencies(s);
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
		while(s->b != NULL)
		{
			gene * t = s->b;
			int gt = (int)fabs(t->x);
			vector<int> v;
			while(true)
			{
				//if(s->x == t->x) break;			// this is a heuristic
				va.push_back(adjacency(PG(s, t)));
				if(t->b == NULL) break;
				if(gc[gt] <= 1) break;
				if(is_fixed_gene(t) == true) break;
				gc[gt]--;
				v.push_back(gt);
				t = t->b;
				gt = (int)fabs(t->x);
			}
			for(int k = 0; k < v.size(); k++) gc[v[k]]++;
			s = s->b;
		}
	}
	return 0;
}

int pbase::build_shared_adjacencies()
{
	vsa.clear();

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
			if(a.e1.g->x == b.e1.g->x && a.e2.g->x == b.e2.g->x && linkable(a.e1.g, b.e1.g) && linkable(a.e2.g, b.e2.g))
			{
				sa = shadj(a.e1.g, a.e2.g, b.e1.g, b.e2.g);
				//if(sa.redundant() == false)
				assert(a.weak_compare(b) >= 2);
				vsa.push_back(sa);
			}

			if(a.e1.g->x + b.e2.g->x == 0 && a.e2.g->x + b.e1.g->x == 0 && linkable(a.e1.g, b.e2.g) && linkable(a.e2.g, b.e1.g))
			{
				sa = shadj(a.e1.g, a.e2.g, b.e2.g, b.e1.g);
				//if(sa.redundant() == false)
				assert(a.weak_compare(b) >= 2);
				vsa.push_back(sa);
			}
		}
	}
	
	return 0;
}

int pbase::remove_redundant_shared_adjacencies()
{
	vector<int> m(vsa.size(), false);
	vector<bool> s(vsa.size(), false);

	int cnt = 0;
	while(true)
	{
		bool f = false;
		for(int i = 0; i < s.size(); i++)
		{
			if(m[i] == true) continue;
			m[i] = true;

			//printf("total = %7d, cnt = %7d\n", (int)s.size(), cnt++);

			bool b = check_redundance(vsa[i], s);

			//string ss = sa.print_string(gi);
			//printf("check %s = %s\n\n", ss.c_str(), b ? "T" : "F");

			if(b == true)
			{
				f = true;
				s[i] = true;
				//break;
			}
		}
		if(f == false) break;
	}

	vector<shadj> v;
	for(int i = 0; i < s.size(); i++)
	{
		if(s[i] == true) continue;
		v.push_back(vsa[i]);
	}

	vsa = v;
	return 0;
}

bool pbase::check_redundance(const shadj & sa, const vector<bool> & s)
{
	vector<gene*> v1;
	vector<gene*> v2;
	v1 = build_gene_list(PG(sa.x1, sa.x2));
	if(sa.direction()) v2 = build_gene_list(PG(sa.y1, sa.y2));
	else v2 = build_gene_list(PG(sa.y2, sa.y1));

	// condition 1: find (x, y), s.t. x1 < x < x2, y1 < y < y2
	// such that (x1, x, y1, y) and (x, x2, y, y2) in s
	map< int, vector<int> > m;
	for(int k = 1; k < v2.size() - 1; k++)
	{
		int f = (int)fabs(v2[k]->x);
		if(m.find(f) == m.end())
		{
			vector<int> v;
			v.push_back(k);
			m.insert(pair< int, vector<int> >(f, v));
		}
		else
		{
			m[f].push_back(k);
		}
	}

	for(int i = 1; i < v1.size() - 1; i++)
	{
		int f = (int)fabs(v1[i]->x);
		if(m.find(f) == m.end()) continue;

		vector<int> & v = m[f];
		for(int k = 0; k < v.size(); k++)
		{
			int j = v[k];
			gene * x = v1[i];
			gene * y = v2[j];
			int p1 = locate_shared_adjacency(sa.x1, x, sa.y1, y);
			int p2 = locate_shared_adjacency(x, sa.x2, y, sa.y2);
			if(p1 == -1) continue;
			if(p2 == -1) continue;
			if(s[p1] == true) continue;
			if(s[p2] == true) continue;
			return true;
		}
	}

	// condition 2: to find x, such that
	// for all (x2, x', y2, y') in s, (x, x', y2, y') in s

	// build pairs that are form shadjs with the two legs of sa
	set<int> s1 = build_pair_contact(sa.x1, sa.y1);
	set<int> s2 = build_pair_contact(sa.x2, sa.y2);
	vector<PG> pv1;
	vector<PG> pv2;
	for(set<int>::iterator it = s1.begin(); it != s1.end(); it++)
	{
		if(s[*it] == true) continue;
		shadj & csa = vsa[*it];
		if(csa.x1 == sa.x1) continue;
		assert(csa.x2 == sa.x1);
		assert(csa.y2 == sa.y1);
		pv1.push_back(PG(csa.x1, csa.y1));
	}
	for(set<int>::iterator it = s2.begin(); it != s2.end(); it++)
	{
		if(s[*it] == true) continue;
		shadj & csa = vsa[*it];
		if(csa.x2 == sa.x2) continue;
		assert(csa.x1 == sa.x2);
		assert(csa.y1 == sa.y2);
		pv2.push_back(PG(csa.x2, csa.y2));
	}

	// four sides
	for(int i = 1; i < v1.size() - 1; i++)
	{
		if(v1[i]->x != sa.x1->x) continue;
		bool f = true;
		for(int k = 0; k < pv1.size(); k++)
		{
			int p = locate_shared_adjacency(pv1[k].first, v1[i], pv1[k].second, sa.y1);
			if(p == -1 || s[p] == true) f = false;
			if(f == false) break;
		}
		if(f == true) return true;
	}
	for(int i = 1; i < v1.size() - 1; i++)
	{
		if(v1[i]->x != sa.x2->x) continue;
		//printf(" test gene %4d:%4d\n", gi[v1[i]], v1[i]->x);

		bool f = true;
		for(int k = 0; k < pv2.size(); k++)
		{
			int p = locate_shared_adjacency(v1[i], pv2[k].first, sa.y2, pv2[k].second);

			//printf(" locate (%4d:%4d) (%4d:%4d) = %3d\n", gi[pv2[k].first], pv2[k].first->x, gi[pv2[k].second], pv2[k].second->x, p);

			if(p == -1 || s[p] == true) f = false;
			if(f == false) break;
		}
		if(f == true) return true;
	}
	for(int i = 1; i < v2.size() - 1; i++)
	{
		if(v2[i]->x != sa.y1->x) continue;
		bool f = true;
		for(int k = 0; k < pv1.size(); k++)
		{
			int p = locate_shared_adjacency(pv1[k].first, sa.x1, pv1[k].second, v2[i]);
			if(p == -1 || s[p] == true) f = false;
			if(f == false) break;
		}
		if(f == true) return true;
	}
	for(int i = 1; i < v2.size() - 1; i++)
	{
		if(v2[i]->x != sa.y2->x) continue;
		bool f = true;
		for(int k = 0; k < pv2.size(); k++)
		{
			int p = locate_shared_adjacency(sa.x2, pv2[k].first, v2[i], pv2[k].second);
			if(p == -1 || s[p] == true) f = false;
			if(f == false) break;
		}
		if(f == true) return true;
	}

	return false;
}

int pbase::build_contact_map()
{
	vct.clear();
	vct.resize(gi.size());
	for(int i = 0; i < vsa.size(); i++)
	{
		shadj & sa = vsa[i];
		vct[gi[sa.x1]].insert(i);
		vct[gi[sa.x2]].insert(i);
		vct[gi[sa.y1]].insert(i);
		vct[gi[sa.y2]].insert(i);
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
			assert(is_fixed_gene(g) == false);
			vsp[gi[g]].insert(i);
			g = g->b;
		}

		p = sort(sa.y1, sa.y2);
		g = p.first->b;
		while(g != p.second)
		{
			assert(gi[g] < vsp.size());
			assert(is_fixed_gene(g) == false);
			vsp[gi[g]].insert(i);
			g = g->b;
		}
	}
	return 0;
}

int pbase::build_candidates()
{
	for(int i = 0; i < vsa.size(); i++) build_candidate(vsa[i]);
	//for(int i = 0; i < cds.size(); i++) cds[i].jp = jumpable(cds[i]);
	return 0;
}

bool pbase::build_candidate(const shadj & sa)
{
	if(sa.adjacent() == false) return false;

	if(is_fixed_gene(sa.x1) == true) return false;
	if(is_fixed_gene(sa.x2) == true) return false;
	if(is_fixed_gene(sa.y1) == true) return false;
	if(is_fixed_gene(sa.y2) == true) return false;

	bool b = sa.direction();
	gene * x = sa.x1;
	gene * y = sa.y1;

	// cannot move left
	if(b == true && x->a != NULL && y->a != NULL && x->a->x == y->a->x)
	{
		if(is_fixed_gene(x->a) == false && is_fixed_gene(y->a) == false) return false;
	}

	if(b == false && x->a != NULL && y->b != NULL && x->a->x + y->b->x == 0)
	{
		if(is_fixed_gene(x->a) == false && is_fixed_gene(y->b) == false) return false;
	}

	vector<gene*> xv;
	vector<gene*> yv;

	while(true)
	{
		xv.push_back(x);
		yv.push_back(y);
		x = x->b;
		if(b) y = y->b;
		else y = y->a;

		if(x == NULL || y == NULL) break;
		if(b == true && x->x != y->x) break;
		if(b == false && x->x + y->x != 0) break;
		if(is_fixed_gene(x) == true) break;
		if(is_fixed_gene(y) == true) break;
	}

	candidate cd(b, xv, yv);
	cds.push_back(cd);

	return true;
}

bool pbase::is_fixed_gene(gene * g)
{
	if(x2y.find(g) != x2y.end()) return true;
	if(y2x.find(g) != y2x.end()) return true;
	return false;
}

bool pbase::is_fixed_pair(gene * x, gene * y)
{
	if(x2y.find(x) == x2y.end()) return false;
	if(x2y[x] != y) return false;
	assert(y2x.find(y) != y2x.end());
	assert(y2x[y] == x);
	return true;
}

bool pbase::linkable(gene * x, gene * y)
{
	int fx = (int)fabs(x->x);
	int fy = (int)fabs(y->x);
	if(fx != fy) return false;
	if(is_fixed_pair(x, y) == true) return true;
	if(is_fixed_gene(x) == true) return false;
	if(is_fixed_gene(y) == true) return false;
	return true;
}

int pbase::add_fixed_pair(gene * x, gene * y)
{
	int fx = (int)fabs(x->x);
	int fy = (int)fabs(y->x);
	assert(fx == fy);

	if(x2y.find(x) == x2y.end() && y2x.find(y) == y2x.end())
	{
		x2y.insert(PG(x, y));
		y2x.insert(PG(y, x));
		vf[fx]++;
	}
	else
	{
		assert(x2y.find(x) != x2y.end());
		assert(y2x.find(y) != y2x.end());
		assert(x2y[x] == y);
		assert(y2x[y] == x);
	}

	return 0;
}

int pbase::fix_singletons()
{
	for(int i = 1; i < gf1.size(); i++)
	{
		if(gf1[i].size() != 1) continue;
		if(gf2[i].size() != 1) continue;
		add_fixed_pair(gf1[i][0], gf2[i][0]);
	}
	return 0;
}

/*
set<int> pbase::set_intersection(const set<int> & x, const set<int> & y)
{
	vector<int> z(x.size() + y.size());
	vector<int>::iterator end = std::set_intersection(x.begin(), x.end(), y.begin(), y.end(), z.begin());
	return set<int>(z.begin(), end);
}

set<int> pbase::set_union(const set<int> & x, const set<int> & y)
{
	vector<int> z(x.size() + y.size());
	vector<int>::iterator end = std::set_union(x.begin(), x.end(), y.begin(), y.end(), z.begin());
	return set<int>(z.begin(), end);
}
*/

set<int> pbase::build_span_intersection(const vector<gene*> & v)
{
	set<int> s;
	if(v.size() == 0) return s;
	int g = gi[v[0]];
	s = vsp[g];
	for(int i = 1; i < v.size(); i++)
	{
		g = gi[v[i]];
		s = set_intersection(s, vsp[g]);
	}
	return s;
}

set<int> pbase::build_span_union(const vector<gene*> & v)
{
	set<int> s;
	for(int i = 0; i < v.size(); i++)
	{
		int g = gi[v[i]];
		s = set_union(s, vsp[g]);
	}
	return s;
}

set<int> pbase::build_contact_intersection(const vector<gene*> & v)
{
	set<int> s;
	if(v.size() == 0) return s;
	int g = gi[v[0]];
	s = vct[g];
	for(int i = 1; i < v.size(); i++)
	{
		int g = gi[v[i]];
		s = set_intersection(s, vct[g]);
	}
	return s;
}


set<int> pbase::build_contact_union(const vector<gene*> & v)
{
	set<int> s;
	for(int i = 0; i < v.size(); i++)
	{
		int g = gi[v[i]];
		s = set_union(s, vct[g]);
	}
	return s;
}

int pbase::locate_shared_adjacency(gene * x1, gene * x2, gene * y1, gene * y2)
{
	vector<gene*> v;
	v.push_back(x1);
	v.push_back(x2);
	v.push_back(y1);
	v.push_back(y2);

	set<int> s = build_contact_intersection(v);

	vector<int> vv;
	for(set<int>::iterator it = s.begin(); it != s.end(); it++)
	{
		shadj & sa = vsa[*it];
		string s = sa.print_string(gi);
		if(sa.x1 == x1 && sa.y1 == y2) continue;
		if(sa.x1 == x2 && sa.y1 == y1) continue;
		vv.push_back(*it);
	}
	assert(vv.size() <= 1);
	if(vv.size() == 0) return -1;
	else if(vv.size() == 1) return vv[0];
	else assert(false);
}

vector<int> pbase::locate_shared_adjacencies(const vector<gene*> & xv, const vector<gene*> & yv)
{
	vector<int> s;
	assert(xv.size() == yv.size());
	if(xv.size() <= 1) return s;

	for(int i = 0; i < xv.size() - 1; i++)
	{
		int p = locate_shared_adjacency(xv[i], xv[i + 1], yv[i], yv[i + 1]);
		assert(p != -1);
		s.push_back(p);
	}
	return s;
}

vector<int> pbase::locate_shared_adjacencies(const candidate & c)
{
	return locate_shared_adjacencies(c.xv, c.yv);
}

bool pbase::innocent(const shadj & sa, const vector<gene*> & xv, const vector<gene*> & yv)
{
	assert(xv.size() == yv.size());

	for(int i = 0; i < xv.size(); i++)
	{
		if(xv[i] == sa.x1 && yv[i] == sa.y1) return true;
		if(xv[i] == sa.x2 && yv[i] == sa.y2) return true;
	}

	return false;
}

set<int> pbase::build_pair_contact(gene * x, gene * y)
{
	set<int> s = set_intersection(vct[gi[x]], vct[gi[y]]);
	set<int> t;
	for(set<int>::iterator it = s.begin(); it != s.end(); it++)
	{
		shadj & sa = vsa[*it];
		if(sa.x1 == x && sa.y1 == y) t.insert(*it);
		if(sa.x2 == x && sa.y2 == y) t.insert(*it);
	}
	return t;
}


