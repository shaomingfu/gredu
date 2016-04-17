#include "pbase.h"

#include <iostream>
#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <cmath>

pbase::pbase(config * _conf, genome * _gm1, genome * _gm2)
	: conf(_conf) 
{
	gm1 = _gm1;
	gm2 = _gm2;
	add_ghost_chrms();

	build_gene_indices();
	build_gene_graph();

	gm1->build_gene_map(gf1);
	gm2->build_gene_map(gf2);
	build_duplicons(gf1, dup1);
	build_duplicons(gf2, dup2);
}

pbase::~pbase()
{}

int pbase::get_index(gene * g)
{
	if(gi1.find(g) != gi1.end()) return gi1[g];
	else if(gi2.find(g) != gi2.end()) return gi2[g] + gi1.size();
	else assert(false);
}

gene * pbase::get_gene(int i)
{
	if(i < ig1.size()) return ig1[i];
	else if(i < gi1.size() + gi2.size()) return ig2[i - gi1.size()];
	else assert(false);
}

int pbase::get_equal_size(int g)
{
	if(out_degree(g, pr) == 0) return 0;
	out_edge_iterator ei1, ei2;
	tie(ei1, ei2) = out_edges(g, pr);
	int h = target(*ei1, pr);
	return out_degree(h, pr);
}

int pbase::get_equal_size(gene * g)
{
	if(gi1.find(g) != gi1.end()) return get_equal_size(gi1[g]);
	else if(gi2.find(g) != gi2.end()) return get_equal_size(gi2[g] + gi1.size());
	else assert(false);
}

int pbase::get_copy_number(gene * g)
{
	int f = (int)fabs(g->x);
	VVG & gf = (gi1.find(g) != gi1.end()) ? gf1 : gf2;
	int n = 0;
	for(int k = 0; k < gf[f].size(); k++)
	{
		gene * h = gf[f][k];
		int x = get_index(h);
		if(out_degree(x, pr) == 0) continue;
		n++;
	}
	return n;
}

int pbase::get_copy_number(int x)
{
	gene * g = get_gene(x);
	return get_copy_number(g);
}

bool pbase::check_all_duplicated(const PG & p)
{
	assert(chrm::order(p.first, p.second));
	vector<gene*> v = build_gene_list(p);
	for(int i = 0; i < v.size(); i++)
	{
		int g = get_index(v[i]);
		if(out_degree(g, pr) <= 0) return false;
		if(get_copy_number(g) <= 1) return false;
	}
	return true;
}


int pbase::add_ghost_chrms()
{
	int n1 = 0;
	int n2 = 0;
	for(int i = 0; i < gm1->chrms.size(); i++) if(gm1->chrms.at(i)->type == LINEAR) n1++;
	for(int i = 0; i < gm2->chrms.size(); i++) if(gm2->chrms.at(i)->type == LINEAR) n2++;

	int n = n1 - n2;
	vector<int> v;
	if(n > 0)
	{
		for(int k = 0; k < n; k++) gm2->add_linear_chrm(v);
	}
	else if(n < 0)
	{
		n = 0 - n;
		for(int k = 0; k < n; k++) gm1->add_linear_chrm(v);
	}
	return 0;
}

int pbase::build_duplicons(const VVG & gf, vector<PG> & dup)
{
	dup.clear();
	for(int f = 1; f < gf.size(); f++)
	{
		if(gf[f].size() <= 1) continue;
		for(int i = 0; i < gf[f].size(); i++)
		{
			int max = 0;
			for(int j = 0; j < gf[f].size(); j++)
			{
				if(j == i) continue;
				gene * x = gf[f][i];
				gene * y = gf[f][j];
				if(x->x == 0) continue;

				int len = 0;
				if(x->x == y->x)
				{
					while(true)
					{
						assert(x != NULL && y != NULL);
						if(x->x == 0 || y->x == 0) break;
						if(x->x != y->x) break;
						len++;
						x = x->b;
						y = y->b;
					}
				}
				else if(x->x + y->x == 0)
				{
					while(true)
					{
						assert(x != NULL && y != NULL);
						if(x->x == 0 || y->x == 0) break;
						if(x->x + y->x != 0) break;
						len++;
						x = x->b;
						y = y->a;
					}
				}
				else assert(false);

				if(len > max) max = len;
			}

			assert(max >= 1);

			gene * x = gf[f][i];
			gene * y = gf[f][i];
			while(max >= 1)
			{
				dup.push_back(PG(x, y));
				y = y->b;
				max--;
			}
		}
	}

	return 0;
}

bool pbase::check_duplicon(const PG & p)
{
	for(int i = 0; i < dup1.size(); i++)
	{
		if(dup1[i] == p) return true;
	}

	for(int i = 0; i < dup2.size(); i++)
	{
		if(dup2[i] == p) return true;
	}

	return false;
}

/*
int pbase::clean_removed_genes()
{
	int n = 0;
	for(int i = 0; i < gi1.size(); i++)
	{
		if(out_degree(i, pr) >= 1) continue;
		n++;
		operation * op = new deletion(ig1[i], ig1[i], true);
		gm1->do_deletion(op);
		delete op;
	}

	for(int i = 0; i < gi2.size(); i++)
	{
		if(out_degree(i + ig1.size(), pr) >= 1) continue;
		n++;
		operation * op = new deletion(ig2[i], ig2[i], true);
		gm2->do_deletion(op);
		delete op;
	}

	return n;
}
*/

int pbase::remove_tandem()
{
	gm1->merge_tandem();
	gm2->merge_tandem();
	return 0;
}

int pbase::build_gene_indices()
{
	gm1->build_gene_indices(gi1, ig1);
	gm2->build_gene_indices(gi2, ig2);
	return 0;
}

int pbase::store_gene_graph(vector<PG> & v)
{
	edge_iterator ei1, ei2;
	for(tie(ei1, ei2) = edges(pr); ei1 != ei2; ei1++)
	{
		int s = source(*ei1, pr);
		int t = target(*ei1, pr);
		gene * x = get_gene(s);
		gene * y = get_gene(t);
		assert(x->x + y->x == 0 || x->x == y->x);
		v.push_back(PG(x, y));
	}
	return 0;
}

int pbase::build_gene_graph(const vector<PG> & v)
{
	pr.clear();

	for(int i = 0; i < ig1.size(); i++) add_vertex(pr);
	for(int i = 0; i < ig2.size(); i++) add_vertex(pr);

	for(int i = 0; i < v.size(); i++)
	{
		gene * x = v[i].first;
		gene * y = v[i].second;
		assert(x->x == y->x || x->x + y->x == 0);
		int s = get_index(x);
		int t = get_index(y);
		add_edge(s, t, pr);
	}
	return 0;
}

int pbase::build_gene_graph()
{
	pr.clear();

	for(int i = 0; i < ig1.size(); i++) add_vertex(pr);
	for(int i = 0; i < ig2.size(); i++) add_vertex(pr);

	vector< vector<gene*> > gf1;
	vector< vector<gene*> > gf2;
	gm1->build_gene_map(gf1);
	gm2->build_gene_map(gf2);

	assert(gf1.size() == gf2.size());
	for(int f = 0; f < gf1.size(); f++)
	{
		for(int i = 0; i < gf1[f].size(); i++)
		{
			for(int j = 0; j < gf2[f].size(); j++)
			{
				gene * g1 = gf1[f][i];
				gene * g2 = gf2[f][j];
				int i1 = gi1[g1];
				int i2 = gi2[g2] + ig1.size();
				add_edge(i1, i2, pr);
			}
		}
	}
	return 0;
}

int pbase::build_edge_map()
{
	ei.clear();
	ie.clear();
	int index = 0;
	edge_iterator ei1, ei2;
	for(tie(ei1, ei2) = edges(pr); ei1 != ei2; ei1++)
	{
		ei.insert(pair<edge_descriptor, int>(*ei1, index));
		ie.insert(pair<int, edge_descriptor>(index, *ei1));
		index++;
	}
	return 0;
}

bool pbase::is_gene_fixed(int x)
{
	if(out_degree(x, pr) != 1) return false;
	out_edge_iterator ei1, ei2;
	tie(ei1, ei2) = out_edges(x, pr);
	int y = target(*ei1, pr);	
	if(out_degree(y, pr) != 1) return false;
	return true;
}

bool pbase::is_gene_pair_fixed(gene * x, gene * y)
{
	if(x == NULL) return false;
	if(y == NULL) return false;
	int s = gi1[x];
	int t = gi2[y] + gi1.size();
	pair<edge_descriptor, bool> p = edge(s, t, pr);
	if(p.second == false) return false;
	if(out_degree(s, pr) >= 2) return false;
	if(out_degree(t, pr) >= 2) return false;
	return true;
}

bool pbase::fix_gene_pair(gene * x, gene * y)
{
	assert(x != NULL);
	assert(y != NULL);
	int a = gi1[x];
	int b = gi2[y] + ig1.size();
	return fix_gene_pair(a, b);
}

bool pbase::fix_gene_pair(int a, int b)
{
	int da = out_degree(a, pr);
	int db = out_degree(b, pr);
	pair<edge_descriptor, bool> p = edge(a, b, pr);
	if(p.second == false) return false;
	if(da == 1 && db == 1) return false;

	clear_vertex(a, pr);
	clear_vertex(b, pr);
	if(a < b) add_edge(a, b, pr);
	else add_edge(b, a, pr);
	return true;
}

int pbase::statistic_degrees()
{
	map<int, int> m;
	for(int i = 0; i < num_vertices(pr); i++)
	{
		int d = out_degree(i, pr);
		if(m.find(d) == m.end()) m.insert(PI(d, 1));
		else m[d]++;
	}

	for(map<int, int>::iterator it = m.begin(); it != m.end(); it++)
	{
		int d = it->first;
		int n = it->second;
		printf("there are %6d vertices of degree %5d\n", n, d);
	}

	return 0;
}

int pbase::statistic_family(const vector<gene*> & v, map<gene*, int> & gi)
{
	for(int i = 0; i < v.size(); i++)
	{
		gene * g = v[i];
		int index = gi[g];
		
		int xaa = (g->a != NULL && g->a->a != NULL) ? g->a->a->x : -1;
		int xa = (g->a != NULL) ? g->a->x : -1;
		int x = g->x;
		int xb = (g->b != NULL) ? g->b->x : -1;
		int xbb = (g->b != NULL && g->b->b != NULL) ? g->b->b->x : -1;
		
		printf("position = %6d, (%7d,%7d,%7d,%7d,%7d)\n", index, xaa, xa, x, xb, xbb);
	}
	return 0;
}

int pbase::statistic_forks()
{
	map<gene*, int>::iterator it;
	for(it = gi1.begin(); it != gi1.end(); it++)
	{
		if(out_degree(it->second, pr) != 2) continue;

		string s;
		string ss = print_triple(it->first, gi1, 0);
		s += "{";
		s += ss;
		s += "}  ";

		out_edge_iterator ei1, ei2;
		int d = -1;
		for(tie(ei1, ei2) = out_edges(it->second, pr); ei1 != ei2; ei1++)
		{
			int t = target(*ei1, pr) - gi1.size();
			if(d == -1) d = out_degree(t + gi1.size(), pr);
			assert(out_degree(t + gi1.size(), pr));
			ss = print_triple(ig2[t], gi2, gi1.size());
			s += "{";
			s += ss;
			s += "}  ";
		}
		if(d != 2) continue;

		printf("%s\n", s.c_str());
	}
	return 0;
}

int pbase::statistic_shared_adjacencies()
{
	edge_iterator ei1, ei2;
	for(tie(ei1, ei2) = edges(pr); ei1 != ei2; ei1++)
	{
		int s = source(*ei1, pr);
		int t = target(*ei1, pr) - ig1.size();

		gene * x = ig1[s];
		gene * y = ig2[t];

		assert(x->x == y->x || x->x + y->x == 0);

		if(x->x == 0) continue;

		if(x->x == y->x && x->b->x == y->b->x)
		{
			int ss = gi1[x->b];
			int tt = gi2[y->b];
			int ds = out_degree(s, pr);
			int dt = out_degree(t + ig1.size(), pr);
			int dss = out_degree(ss, pr);
			int dtt = out_degree(tt + ig1.size(), pr);
			int d = ds + dt + dss + dtt;
			if(d == 4) continue;
			printf("shared adjacency, gene = (%5d,%5d), (%5d,%5d), position = (%5d,%5d), (%5d,%5d), ", x->x, x->b->x, y->x, y->b->x, s, ss, t, tt);
			printf("degree = (%2d,%2d), (%2d,%2d), total = %3d\n", ds, dss, dt, dtt, d);
		}
		else if(x->x + y->x == 0 && x->b->x + y->a->x == 0)
		{
			int ss = gi1[x->b];
			int tt = gi2[y->a];
			int ds = out_degree(s, pr);
			int dt = out_degree(t + ig1.size(), pr);
			int dss = out_degree(ss, pr);
			int dtt = out_degree(tt + ig1.size(), pr);
			int d = ds + dt + dss + dtt;
			if(d == 4) continue;
			printf("shared adjacency, gene = (%5d,%5d), (%5d,%5d), position = (%5d,%5d), (%5d,%5d), ", x->x, x->b->x, y->x, y->a->x, s, ss, t, tt);
			printf("degree = (%2d,%2d), (%2d,%2d), total = %3d\n", ds, dss, dt, dtt, d);
		}
	}

	return 0;
}

int pbase::statistic_fixed_gene_pairs()
{
	int s1 = 0;
	int s2 = 0;
	for(int i = 0; i < gi1.size(); i++)
	{
		if(is_gene_fixed(i)) s1++;
	}

	for(int i = 0; i < gi2.size(); i++)
	{
		if(is_gene_fixed(i + gi1.size())) s2++;
	}

	assert(s1 == s2);
	printf("total %6d fixed gene pairs\n", s1);
	return s1;
}

int pbase::build_fixed_gene_pairs()
{
	px2y.clear();
	py2x.clear();
	for(int x = 0; x < gi1.size(); x++)
	{
		if(out_degree(x, pr) != 1) continue;
		out_edge_iterator ei1, ei2;
		tie(ei1, ei2) = out_edges(x, pr);
		int y = target(*ei1, pr);	
		if(out_degree(y, pr) != 1) continue;

		gene * g = ig1[x];
		gene * h = ig2[y - gi1.size()];
		px2y.insert(PG(g, h));
		py2x.insert(PG(h, g));
	}
	return 0;
}


int pbase::print_gene_graph()
{
	edge_iterator ei1, ei2;
	for(tie(ei1, ei2) = edges(pr); ei1 != ei2; ei1++)
	{
		int s = source(*ei1, pr);
		int t = target(*ei1, pr);
		int ds = out_degree(s, pr);
		int dt = out_degree(t, pr);
		gene * x = ig1[s];
		gene * y = ig2[t - ig1.size()];
		if(x->x == 0) continue;
		printf("%20s %20s %3d %3d\n", x->s.c_str(), y->s.c_str(), ds, dt);
	}
	return 0;
}

string pbase::print_triple(gene * x, map<gene*, int> & gi, int n)
{
	gene * a = x->a;
	gene * b = x->b;
	int dx = out_degree(gi[x] + n, pr);
	int da = out_degree(gi[a] + n, pr);
	int db = out_degree(gi[b] + n, pr);
	char buf[10240];
	sprintf(buf, "(%6d,%6d,%6d), (%2d,%2d,%2d)", a->x, x->x, b->x, da, dx, db);
	return string(buf);
}

