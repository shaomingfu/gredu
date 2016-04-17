#include "asolver.h"
#include <iostream>
#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <cmath>

asolver::asolver(config * _conf, genome * _gm1, genome * _gm2)
	: pbase(_conf, _gm1, _gm2)
{
	statistic_fixed_gene_pairs();

	build_complements();
	build_adjacencies();
	build_adjacency_graph();

	printf("total %6d edges in gene graph\n", (int)num_edges(pr));
	while(true)
	{
		resolve_adjacencies();
		int n1 = identify_shared_adjacencies();
		int n2 = identify_duplicated_circulars(dup1);
		int n3 = identify_duplicated_circulars(dup2);
		int n4 = identify_ltypes();
		int n5 = identify_vtypes();
		int n6 = identify_8types();
		printf("identify (%4d, %4d + %4d, %4d, %4d, %4d) substructures, total %6d edges in gene graph\n", n1, n2, n3, n4, n5, n6, (int)num_edges(pr));
		if(n1 + n2 + n3 + n4 + n5 + n6 == 0) break;
	}

	statistic_fixed_gene_pairs();
	//build_fixed_gene_pairs();

	/*
	printf("genome1:\n");
	print_genome(gm1);
	printf("genome2:\n");
	print_genome(gm2);
	*/
}

asolver::~asolver()
{}

int asolver::build_complements()
{
	cpl.assign(num_vertices(pr) * 2, -1);
	for(int i = 0; i < num_vertices(pr); i++)
	{
		int x = i * 2 + 0;
		int y = i * 2 + 1;
		cpl[x] = y;
		cpl[y] = x;
	}
	return 0;
}

int asolver::build_adjacency_graph()
{
	ar.clear();
	for(int i = 0; i < num_vertices(pr); i++)
	{
		add_vertex(ar);
		add_vertex(ar);
	}

	edge_iterator ei1, ei2;
	for(tie(ei1, ei2) = edges(pr); ei1 != ei2; ei1++)
	{
		int s = source(*ei1, pr);
		int t = target(*ei1, pr);
		add_edge(s * 2 + 0, t * 2 + 0, ar);
		add_edge(s * 2 + 1, t * 2 + 1, ar);
	}

	return 0;
}

int asolver::build_adjacencies()
{
	adj.assign(num_vertices(pr) * 2, -1);
	map<gene*, int>::iterator it;
	for(it = gi1.begin(); it != gi1.end(); it++)
	{
		gene * g = it->first;
		gene * h = g->b;
		if(h == NULL) continue;
		PI p = make_adjacency(g, h);
		adj[p.first] = p.second;
		adj[p.second] = p.first;
	}

	for(it = gi2.begin(); it != gi2.end(); it++)
	{
		gene * g = it->first;
		gene * h = g->b;
		if(h == NULL) continue;
		PI p = make_adjacency(g, h);
		adj[p.first] = p.second;
		adj[p.second] = p.first;
	}

	for(int i = 0; i < adj.size(); i++)
	{
		int j = adj[i];
		if(j == -1) continue;
		assert(adj[i] == j);
		assert(adj[j] == i);
	}
	return 0;
}

int asolver::identify_duplicated_circulars(const vector<PG> & dup)
{
	int cnt = 0;
	for(int i = 0; i < dup.size(); i++)
	{
		bool b = check_duplicated_circular(dup[i]);
		if(b == false) continue;
		remove_segment(dup[i]);
		cnt++;
	}
	resolve_removed_genes();
	return cnt++;
}

bool asolver::check_duplicated_circular(const PG & p)
{
	if(check_availability(p) == false) return false;
	if(adjacent(p.second, p.first) == false) return false;
	if(check_all_duplicated(p) == false) return false;
	return true;
}

bool asolver::check_availability(const PG & p)
{
	assert(p.first != NULL);
	assert(p.second != NULL);
	gene * g = p.first;
	while(true)
	{
		int x = get_index(g);
		if(out_degree(x, pr) == 0) return false;
		if(g == p.second) return true;
		if(g->b == NULL) return false;
		if(adjacent(g, g->b) == false) return false;
		g = g->b;
	}
	assert(false);
}

int asolver::remove_segment(const PG & p)
{
	assert(check_availability(p));
	assert(check_all_duplicated(p));
	vector<gene*> v = build_gene_list(p);
	for(int i = 0; i < v.size(); i++)
	{
		int x = get_index(v[i]);
		clear_vertex(x, pr);
		clear_vertex(x * 2 + 0, ar);
		clear_vertex(x * 2 + 1, ar);
	}
	return 0;
}

int asolver::resolve_adjacencies()
{
	// TODO
	int cnt = 0;
	while(true)
	{
		bool b = false;
		for(int i = 0; i < adj.size(); i++)
		{
			int j = adj[i];
			if(i > j) continue;
			bool f = resolve_adjacency(i, j);
			if(f == true) b = true;
			if(f == true) cnt++;
		}
		if(b == false) break;
	}
	return cnt;
}

bool asolver::resolve_adjacency(int x, int y)
{
	assert(x < y);
	assert(adj[x] == y);
	assert(adj[y] == x);
	if(out_degree(x, ar) != 1) return false;
	if(out_degree(y, ar) != 1) return false;
	out_edge_iterator ei1, ei2;
	tie(ei1, ei2) = out_edges(x, ar);
	int s = target(*ei1, ar);
	tie(ei1, ei2) = out_edges(y, ar);
	int t = target(*ei1, ar);
	if(adj[s] == t) return false;
	if(out_degree(s, ar) != 1) return false;
	if(out_degree(t, ar) != 1) return false;
	assert(adj[t] != s);
	int a = adj[s];
	int b = adj[t];
	assert(adj[a] == s);
	assert(adj[b] == t);
	adj[a] = b;
	adj[b] = a;
	adj[s] = t;
	adj[t] = s;
	return true;
}

int asolver::resolve_removed_genes()
{
	// TODO
	for(int i = 0; i < num_vertices(pr); i++)
	{
		if(out_degree(i, pr) >= 1) continue;
		if(get_gene(i) == NULL) continue;
		resolve_removed_gene(i);
	}
	return 0;
}

int asolver::resolve_removed_gene(int i)
{
	assert(out_degree(i, pr) == 0);
	int x = i * 2 + 0;
	int y = i * 2 + 1;
	int dx = out_degree(x, ar);
	int dy = out_degree(y, ar);
	gene * g = get_gene(i);	
	//printf("remove gene %5d = %5d, (%5d, %5d), degree = (%2d, %2d)\n", i, g->x, x, y, dx, dy);
	assert(dx == 0);
	assert(dy == 0);
	int s = adj[x];
	int t = adj[y];
	adj[x] = y;
	adj[y] = x;
	adj[s] = t;
	adj[t] = s;
	return 0;
}

int asolver::identify_shared_adjacencies()
{
	int cnt = 0;
	edge_iterator ei1, ei2;
	for(tie(ei1, ei2) = edges(ar); ei1 != ei2; ei1++)
	{
		int x = source(*ei1, ar);
		int y = target(*ei1, ar);
		if(x > y) continue;
		bool f = identify_shared_adjacency(x, y);
		if(f == true) cnt++;
	}
	resolve_removed_genes();
	return cnt;
}

bool asolver::identify_shared_adjacency(int x, int y)
{
	if(out_degree(x, ar) != 1) return false;
	if(out_degree(y, ar) != 1) return false;
	assert(edge(x, y, ar).second == true);
	int s = adj[x];
	int t = adj[y];

	if(s == -1) return false;
	assert(t != -1);

	if(edge(s, t, ar).second == false) return false;
	int a = cpl[s];
	int b = cpl[t];
	assert(edge(a, b, ar).second == true);
	int g = s / 2;
	int h = t / 2;
	assert(g == a / 2);
	assert(h == b / 2);
	bool b1 = fix_extremity_pair(s, t);
	bool b2 = fix_extremity_pair(a, b);
	bool b3 = fix_gene_pair(g, h);
	assert(b1 == b2);
	assert(b1 == b3);
	return b1;
}

int asolver::identify_ltypes()
{
	int cnt = 0;
	edge_iterator ei1, ei2;
	for(tie(ei1, ei2) = edges(ar); ei1 != ei2; ei1++)
	{
		int x = source(*ei1, ar);
		int y = target(*ei1, ar);
		assert(x != y);
		if(x > y) continue;
		bool f = identify_ltype(x, y);
		if(f == true) cnt++;
	}
	resolve_removed_genes();
	return cnt;
}

bool asolver::identify_ltype(int x, int y)
{
	if(out_degree(x, ar) != 1) return false;
	if(out_degree(y, ar) != 1) return false;
	assert(edge(x, y, ar).second == true);
	int s = adj[x];
	int t = adj[y];

	if(s == -1) return false;
	assert(t != -1);

	int cs = get_equal_size(s / 2);
	int ct = get_equal_size(t / 2);

	if(cs != 1 && ct != 1) return false;
	if(cs == 1 && ct != 1) return identify_ltype(y, x);
	
	assert(ct == 1);

	bool b = false;
	PG p(NULL, NULL);
	out_edge_iterator ei1, ei2;
	int w = -1;
	for(tie(ei1, ei2) = out_edges(t, ar); ei1 != ei2; ei1++)
	{
		w = target(*ei1, ar);
		if(adj[w] == x) continue;
		if(cpl[w] == x) continue;
		if(dcj_direction(x, w) == false) continue;
		PG q = make_adjacency(x, w);
		if(check_availability(q) == false) continue;
		p = PG(q.first->b, q.second->a);
		if(check_all_duplicated(p) == false) continue;
		if(check_duplicon(p) == false) continue;
		b = true;
		break;
	}

	if(b == false) return false;

	remove_segment(p);

	int tt = cpl[t];
	int ww = cpl[w];
	assert(edge(ww, tt, ar).second == true);
	bool b1 = fix_extremity_pair(w, t);
	bool b2 = fix_extremity_pair(ww, tt);
	bool b3 = fix_gene_pair(w / 2, t / 2);
	assert(b1 == b2);
	assert(b1 == b3);

	return true;
}

int asolver::identify_vtypes()
{
	int cnt = 0;
	for(int i = 0; i < num_vertices(ar); i++)
	{
		int j = adj[i];
		if(j == -1) continue;
		if(i > j) continue;
		if(cpl[i] == j) continue;
		bool b = identify_vtype(i, j);
		if(b == true) cnt++;
	}
	resolve_removed_genes();
	return cnt;
}

bool asolver::identify_vtype(int x, int y)
{
	assert(adj[x] == y);
	assert(adj[y] == x);
	assert(cpl[x] != y);
	assert(cpl[y] != x);

	if(get_copy_number(x / 2) >= 2) return false;
	if(get_copy_number(y / 2) >= 2) return false;

	int dx = out_degree(x, ar);
	int dy = out_degree(y, ar);
	if(dx == 1 && dy == 1) return false;

	if(shared(x, y) == false) return false;

	bool bxy = dominate(x, y);
	bool byx = dominate(y, x);

	if(bxy == byx) return false;
	if(byx == true) return identify_vtype(y, x);
	assert(bxy == true);
	assert(byx == false);
	assert(dy > dx);

	out_edge_iterator ei1, ei2;
	vector<int> v;
	for(tie(ei1, ei2) = out_edges(y, ar); ei1 != ei2; ei1++)
	{
		int s = target(*ei1, ar);
		int t = adj[s];
		if(edge(t, x, ar).second) continue;
		v.push_back(s);
	}

	assert(v.size() >= 1);

	for(int i = 0; i < v.size(); i++)
	{
		int s = v[i];
		assert(edge(s, y, ar).second);
		remove_edge(s, y, ar);
		assert(edge(cpl[s], cpl[y], ar).second);
		remove_edge(cpl[s], cpl[y], ar);
		assert(edge(s / 2, y / 2, pr).second);
		remove_edge(s / 2, y / 2, pr);
	}

	return true;
}

int asolver::identify_8types()
{
	int cnt = 0;
	for(int i = 0; i < num_vertices(ar); i++)
	{
		int j = adj[i];
		if(j == -1) continue;
		if(i > j) continue;
		if(cpl[i] == j) continue;
		bool b = identify_8type(i, j);
		if(b == true) cnt++;
	}
	resolve_removed_genes();
	return cnt;
}

bool asolver::identify_8type(int x, int y)
{
	assert(adj[x] == y);
	assert(adj[y] == x);
	assert(cpl[x] != y);
	assert(cpl[y] != x);

	if(get_copy_number(x / 2) >= 2) return false;
	if(get_copy_number(y / 2) >= 2) return false;

	int dx = out_degree(x, ar);
	int dy = out_degree(y, ar);
	if(dx == 1 && dy == 1) return false;
	
	out_edge_iterator ei1, ei2;
	out_edge_iterator ei3, ei4;
	vector<PI> vp;
	vector< vector<PI> > vc;
	for(tie(ei1, ei2) = out_edges(x, ar); ei1 != ei2; ei1++)
	{
		int s = target(*ei1, ar);
		int a = adj[s];
		if(a == -1) continue;
		for(tie(ei3, ei4) = out_edges(y, ar); ei3 != ei4; ei3++)
		{
			int t = target(*ei3, ar);	
			int b = adj[t];
			if(b == -1) continue;
			vector<PI> v = find_possible_adjacencies(a, b);
			vp.push_back(PI(s, t));
			vc.push_back(v);
		}
	}

	int n = 0;
	int k = -1;
	for(int i = 0; i < vp.size(); i++) 
	{
		if(vc[i].size() == 0) continue;
		n++;
		k = i;
	}
	
	if(n != 1) return false;
	if(vc[k].size() != 1) return false;

	//printf("8type (%6d, %6d): %3d possible adjacencies\n", x, y, m);

	int s = vp[k].first;
	int t = vp[k].second;
	int a = adj[s];
	int b = adj[t];
	int p = vc[k][0].first;
	int q = vc[k][0].second;
	assert(adj[p] == q);
	assert(adj[q] == p);

	if(b == cpl[a] || b == cpl[s]) return false;
	if(t == cpl[a] || t == cpl[s]) return false;
	if(p == cpl[x] || p == cpl[y]) return false;
	if(q == cpl[x] || q == cpl[y]) return false;

	/*
	printf("(s, x) = (%5d,%5d), (t, y) = (%5d,%5d), (a, p) = (%5d,%5d), (b, q) = (%5d,%5d)\n", s, x, t, y, a, p, b, q);
	printf("(s, x) = (%5d,%5d), (t, y) = (%5d,%5d), (a, p) = (%5d,%5d), (b, q) = (%5d,%5d)\n", cpl[s], cpl[x], cpl[t], cpl[y], cpl[a], cpl[p], cpl[b], cpl[q]);
	printf("\n");
	*/

	fix_gene_pair(s / 2, x / 2);
	fix_gene_pair(t / 2, y / 2);
	fix_gene_pair(a / 2, p / 2);
	fix_gene_pair(b / 2, q / 2);

	fix_extremity_pair(s, x);
	fix_extremity_pair(t, y);
	fix_extremity_pair(a, p);
	fix_extremity_pair(b, q);
	fix_extremity_pair(cpl[s], cpl[x]);
	fix_extremity_pair(cpl[t], cpl[y]);
	fix_extremity_pair(cpl[a], cpl[p]);
	fix_extremity_pair(cpl[b], cpl[q]);

	return true;
}

vector<PI> asolver::find_possible_adjacencies(int x, int y)
{
	out_edge_iterator ei1, ei2;
	out_edge_iterator ei3, ei4;
	vector<PI> v;
	for(tie(ei1, ei2) = out_edges(x, ar); ei1 != ei2; ei1++)
	{
		int s = target(*ei1, ar);
		for(tie(ei3, ei4) = out_edges(y, ar); ei3 != ei4; ei3++)
		{
			int t = target(*ei3, ar);
			if(adj[s] != t) continue;
			assert(adj[t] == s);
			v.push_back(PI(s, t));
		}
	}
	return v;
}

bool asolver::fix_extremity_pair(int a, int b)
{
	assert(edge(a, b, ar).second);
	int da = out_degree(a, ar);
	int db = out_degree(b, ar);
	if(da == 1 && db == 1) return false;
	clear_vertex(a, ar);
	clear_vertex(b, ar);
	if(a < b) add_edge(a, b, ar);
	else add_edge(b, a, ar);
	return true;
}

int asolver::print_adjacency_graph()
{
	for(int i = 0; i < num_vertices(ar); i++)
	{
		printf("extremity %5d, complement = %5d, adjacenct extremity = %5d, degree = %5d\n", i, cpl[i], adj[i], (int)out_degree(i, ar));
	}
	return 0;
}

int asolver::print_genome(genome * gm)
{
	for(int i = 0; i < gm->chrms.size(); i++)
	{
		chrm * ch = gm->chrms.at(i);
		if(ch->type == LINEAR)
		{
			gene * q = ch->p->b;
			while(q != NULL)
			{
				PI p = make_extremities(q);
				int g = get_index(q);
				printf("[L]chr%2d, gene%5d=%7d, (%5d,%5d), equal = %2d, copy = %2d, degree = %2d\n",
						i + 1, g, q->x, p.first, p.second, get_equal_size(g), get_copy_number(g), int(out_degree(g, pr)));
				q = q->b;
			}
		}
		else
		{
			assert(ch->p != NULL);
			gene * q = ch->p;
			while(true)
			{
				PI p = make_extremities(q);
				int g = get_index(q);
				printf("[C]chr%2d, gene%5d=%7d, (%5d,%5d), equal = %2d, copy = %2d, degree = %2d\n",
						i + 1, g, q->x, p.first, p.second, get_equal_size(g), get_copy_number(g), int(out_degree(g, pr)));
				q = q->b;
				if(q == ch->p) break;
			}
		}
	}
	return 0;
}

PG asolver::make_adjacency(int a, int b)
{
	assert(dcj_direction(a, b));
	gene * x = get_gene(a / 2);
	gene * y = get_gene(b / 2);
	PI s = make_extremities(x);
	PI t = make_extremities(y);

	//printf("(a, b) = (%5d, %5d), x = %6d (%5d, %5d), y = %6d (%5d, %5d)\n", a, b, x->x, s.first, s.second, y->x, t.first, t.second);
	if(a == s.second && b == t.first) return PG(x, y);
	else if(a == s.first && b == t.second) return PG(y, x);
	else assert(false);
}

PI asolver::make_adjacency(gene * a, gene * b)
{
	PI x = make_extremities(a);
	PI y = make_extremities(b);
	return PI(x.second, y.first);
}

PI asolver::make_extremities(gene * g)
{
	int x = get_index(g);
	int a = 2 * x + 0;			// tail
	int b = 2 * x + 1;			// head
	if(g->x == 0 && g->a == NULL) return PI(a, b);
	else if(g->x == 0 && g->b == NULL) return PI(b, a);
	else if(g->x > 0) return PI(a, b);
	else if(g->x < 0) return PI(b, a);
	else assert(false);
}

bool asolver::adjacent(gene * x, gene * y)
{
	PI p = make_adjacency(x, y);
	int s = p.first;
	int t = p.second;
	if(adj[s] != t) return false;
	assert(adj[t] == s);
	return true;
}

bool asolver::dcj_direction(int x, int y)
{
	PI a = make_extremities(get_gene(x / 2));
	PI b = make_extremities(get_gene(y / 2));
	if(x == a.second && y == b.first) return true;
	if(x == a.first && y == b.second) return true;
	return false;
}

bool asolver::shared(int x, int y)
{	
	assert(adj[x] == y);
	assert(adj[y] == x);
	vector<PI> v = find_possible_adjacencies(x, y);
	if(v.size() >= 1) return true;
	else return false;
}

bool asolver::dominate(int x, int y)
{	
	assert(adj[x] == y);
	assert(adj[y] == x);
	out_edge_iterator ei1, ei2;
	for(tie(ei1, ei2) = out_edges(x, ar); ei1 != ei2; ei1++)
	{
		int s = target(*ei1, ar);
		int t = adj[s];
		if(edge(t, y, ar).second == false) return false;
	}
	return true;
}

