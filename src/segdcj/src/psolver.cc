#include "psolver.h"
#include <boost/graph/connected_components.hpp>
#include <iostream>
#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <cmath>

psolver::psolver(config * _conf, genome * _gm1, genome * _gm2)
	: bsolver(_conf, _gm1, _gm2)
{
	// rebuild all the data structures
	/*
	vector<PG> v;
	store_gene_graph(v);
	build_gene_indices();
	build_gene_graph(v);
	build_edge_map();

	build_complements();
	build_adjacencies();
	build_adjacency_graph();

	gm1->build_gene_map(gf1);
	gm2->build_gene_map(gf2);
	pbase::build_duplicons(gf1, dup1);
	pbase::build_duplicons(gf2, dup2);
	*/

	build_edge_map();
	build_gene_list();
	build_gene_families();

	transform_duplicons();

	//build_inference_pairs();
}

psolver::~psolver()
{}

int psolver::build_gene_list()
{
	ge.clear();
	for(int i = 0; i < num_vertices(pr); i++)
	{
		ge.push_back(PI(i * 2 + 0, i * 2 + 1));
	}
	return 0;
}

int psolver::build_gene_families()
{
	gf.assign(num_vertices(pr), -1);
	for(int i = 0; i < gf.size(); i++)
	{
		gene * g = get_gene(i);
		int f = (int)fabs(g->x);
		if(i >= gi1.size()) f = 0 - f;
		gf[i] = f;
	}
	return 0;
}

int psolver::transform_duplicons()
{
	vdup.clear();
	for(int i = 0; i < dup1.size(); i++)
	{
		gene * x = dup1[i].first;
		gene * y = dup1[i].second;
		vector<int> v;
		while(true)
		{
			int g = get_index(x);
			v.push_back(g);
			if(x == y) break;
			x = x->b;
		}
		vdup.push_back(v);
	}

	for(int i = 0; i < dup2.size(); i++)
	{
		gene * x = dup2[i].first;
		gene * y = dup2[i].second;
		vector<int> v;
		while(true)
		{
			int g = get_index(x);
			v.push_back(g);
			if(x == y) break;
			x = x->b;
		}
		vdup.push_back(v);
	}

	return 0;
}

int psolver::build_inference_pairs()
{
	edge_iterator ei1, ei2;
	for(tie(ei1, ei2) = edges(pr); ei1 != ei2; ei1++)
	{
		int s = source(*ei1, pr);
		int t = target(*ei1, pr);

		gene * x = ig1[s];
		gene * y = ig2[t - ig1.size()];

		assert(x->x == y->x || x->x + y->x == 0);

		if(x->x == 0) continue;

		int f = (int)fabs(x->x);

		if(x->x == y->x && x->b->x == y->b->x) add_inference_pair(s, t, x->b, y->b);
		if(x->x + y->x == 0 && x->b->x + y->a->x == 0) add_inference_pair(s, t, x->b, y->a);
	}

	//printf("total %6d inference pairs\n", (int)vip.size());

	return 0;
}

bool psolver::add_inference_pair(int s, int t, gene * x, gene * y)
{
	assert(x != NULL);
	assert(y != NULL);
	if(x->x == 0) return false;
	if(y->x == 0) return false;
	assert(x->x == y->x || x->x + y->x == 0);
	int a = gi1[x];
	int b = gi2[y] + ig1.size();
	pair<edge_descriptor, bool> p1 = edge(s, t, pr);
	pair<edge_descriptor, bool> p2 = edge(a, b, pr);
	assert(p1.second == true);
	if(p2.second == false) return false;
	assert(edge(s, t, pr).second == true);
	int da = out_degree(a, pr);
	int db = out_degree(b, pr);
	int ds = out_degree(s, pr);
	int dt = out_degree(t, pr);
	if(da==1 && db==1 && ds==1 && dt==1) return false;
	int f = (int)fabs(x->x);
	assert(f == (int)fabs(y->x));
	vip.push_back(PI(ei[p1.first], ei[p2.first]));
	return true;
}

map<gene*, gene*> psolver::transform_mapping(const map<int, int> & m)
{
	map<gene*, gene*> x;
	for(map<int, int>::const_iterator it = m.begin(); it != m.end(); it++)
	{
		int s = it->first;
		int t = it->second;
		gene * g = get_gene(s);
		gene * h = get_gene(t);
		x.insert(PG(g, h));
	}
	return x;
}

int psolver::check_families()
{
	for(int i = 0; i < gf1.size(); i++)
	{
		int n1 = 0;
		for(int j = 0; j < gf1[i].size(); j++)
		{
			gene * g = gf1[i][j];
			int x = gi1[g];
			int d = out_degree(x, pr);
			if(d >= 1) n1++;
		}
		int n2 = 0;
		for(int j = 0; j < gf2[i].size(); j++)
		{
			gene * g = gf2[i][j];
			int x = gi2[g] + gi1.size();
			int d = out_degree(x, pr);
			if(d >= 1) n2++;
		}
		assert(n1 >= 1 && n2 >= 1);
		printf("(family %5d, size = (%4d, %4d)\n", i, n1, n2);
	}
	return 0;
}
