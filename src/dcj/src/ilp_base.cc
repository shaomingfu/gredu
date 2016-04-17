#include "ilp_base.h"
#include "mygraph.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <boost/graph/connected_components.hpp>


using namespace std;

ilp_base::ilp_base(const ugraph & _gr, const vector<PII> & _gene_list, const vector<int> & _partners)
	: gr(_gr), gene_list(_gene_list), partners(_partners)
{
	env = new GRBEnv();
	model = new GRBModel(*env);
    //model->getEnv().set(GRB_IntParam_OutputFlag, 0);

	nodes = 10000000;
	gap = 1e-4;
	cuts = -1;
	prepasses = -1;
	heuristics = 0.05;
	mipfocus = 0;
	gomorypasses = -1;
	scaleflag = 1;
	numericfocus = 0;

	//verify_bipartite_graph();
}

ilp_base::~ilp_base()
{
	delete model;
	delete env;
}

int ilp_base::verify_bipartite_graph()
{
	assert(num_vertices(gr) % 4 == 0);
	int n = num_vertices(gr) / 2;
	
	out_edge_iterator ei1, ei2;
	for(int i = 0; i < n; i++)
	{
		for(tie(ei1, ei2) = out_edges(i, gr); ei1 != ei2; ei1++)
		{
			assert(source(*ei1, gr) == i);
			assert(target(*ei1, gr) >= n);
		}
	}

	for(int i = n; i < 2 * n; i++)
	{
		for(tie(ei1, ei2) = out_edges(i, gr); ei1 != ei2; ei1++)
		{
			assert(source(*ei1, gr) == i);
			assert(target(*ei1, gr) <  n);
		}
	}

	for(int i = 0; i < gene_list.size(); i++)
	{
		if(gene_list.at(i).second == 0) assert(gene_list.at(i).first.second == gene_list.at(i).first.first);
		else if(gene_list.at(i).first.first < n) assert(gene_list.at(i).first.second < n);
		else if(gene_list.at(i).first.first >= n) assert(gene_list.at(i).first.second >= n);
		else assert(false);
		assert(gene_list.at(i).second >= 0);
	}

	assert(partners.size() == 2 * n);
	for(int i = 0; i < partners.size(); i++)
	{
		if(i < n) assert(partners[i] < n);
		if(i >= n) assert(partners[i] >= n);
	}

	return 0;
}

int ilp_base::build_neighbor_indices()
{
	int n = num_vertices(gr) / 2;
	neighbors.clear();
	nb_indices.clear();
	nb_indices.assign(2 * n, -1);
	for(int i = 0; i < 2 * n; i++)
	{
		if(nb_indices[i] != -1) continue;
		if(i < partners[i]) neighbors.push_back(PI(i, partners[i]));
		else neighbors.push_back(PI(partners[i], i));
		nb_indices[i] = neighbors.size() - 1;
		nb_indices[partners[i]] = neighbors.size() - 1;
	}
	assert(neighbors.size() == n);
	return 0;
}

int ilp_base::build_gene_indices()
{
	gn_indices.assign(num_vertices(gr), -1);
	for(int i = 0; i < gene_list.size(); i++)
	{
		gn_indices[gene_list[i].first.first] = i;
		gn_indices[gene_list[i].first.second] = i;
	}
	return 0;
}

int ilp_base::build_complements()
{
	complements.assign(num_vertices(gr), -1);
	for(int i = 0; i < gene_list.size(); i++)
	{
		complements[gene_list[i].first.first] = gene_list[i].first.second;
		complements[gene_list[i].first.second] = gene_list[i].first.first;
	}
	return 0;
}

int ilp_base::build_gene_pairs()
{
	int n = num_vertices(gr) / 2;
	gene_pairs.clear();
	pair<edge_descriptor, bool> p1, p2;
	for(int i = 0; i < gene_list.size(); i++)
	{
		PII x = gene_list.at(i);
		int x1 = x.first.first;
		int x2 = x.first.second;
		for(int j = i + 1; j < gene_list.size(); j++)
		{
			PII y = gene_list.at(j);
			int y1 = y.first.first;
			int y2 = y.first.second;

			if(x1 < n && y1 < n) continue;
			if(x1 >= n && y1 >= n) continue;
			if(x.second != y.second) continue;

			p1 = edge(x1, y1, gr);
			p2 = edge(x2, y2, gr);

			//printf("x1 = %4d, dx = %4d, y1 = %4d, dy = %4d, gene = %4d\n", x1, out_degree(x1, gr), y1, out_degree(y1, gr), gene_list[i].second);
			//printf("x2 = %4d, dx = %4d, y2 = %4d, dy = %4d, gene = %4d\n", x2, out_degree(x2, gr), y2, out_degree(y2, gr), gene_list[i].second);

			assert(p1.second == p2.second);
			if(p1.second == false) continue;

			/*
			printf("x1 = %4d, dx = %4d, y1 = %4d, dy = %4d, gene = %4d\n", x1, out_degree(x1, gr), y1, out_degree(y1, gr), gene_list[i].second);
			assert(out_degree(x1, gr) == out_degree(y1, gr));
			assert(out_degree(x2, gr) == out_degree(y2, gr));
			*/

			//printf("x:(%4d,%4d,%4d), y:(%4d,%4d,%4d)\n", x1, x2, x.second, y1, y2, y.second);

			assert(out_degree(x1, gr) == out_degree(x2, gr));
			assert(out_degree(y1, gr) == out_degree(y2, gr));
			if(out_degree(x1, gr) == 1 || out_degree(y1, gr) == 1)
			{
				pairs.insert(PI(x1, y1));
				pairs.insert(PI(x2, y2));
				//assert(pairs.find(PI(x1, y1)) != pairs.end());
				//assert(pairs.find(PI(x2, y2)) != pairs.end());
				continue;
			}

			if(x1 < n && y1 >= n) gene_pairs.push_back(PI(i, j));
			if(x1 >= n && y1 < n) gene_pairs.push_back(PI(j, i));
		}
	}
	return 0;
}

int ilp_base::build_edge_map()
{
	edge_iterator ei1, ei2;
	for(tie(ei1, ei2) = edges(gr); ei1 != ei2; ++ei1)
	{
		eim[*ei1] = -1;
	}

	for(int k = 0; k < gene_pairs.size(); k++)
	{
		PII g1 = gene_list[gene_pairs[k].first];
		PII g2 = gene_list[gene_pairs[k].second];

		//printf("gene pair %4d: (%4d, %4d), (%4d, %4d)\n", k, g1.first.first, g1.first.second, g2.first.first, g2.first.second);

		pair<edge_descriptor, bool> peb1 = edge(g1.first.first, g2.first.first, gr);
		pair<edge_descriptor, bool> peb2 = edge(g1.first.second, g2.first.second, gr);

		assert(peb1.second == true);
		assert(peb2.second == true);

		assert(eim.find(peb1.first) != eim.end());
		assert(eim.find(peb2.first) != eim.end());

		eim[peb1.first] = k;
		eim[peb2.first] = k;
	}

	return 0;
}

int ilp_base::set_timelimit(double timelimit)
{
	model->getEnv().set(GRB_DoubleParam_TimeLimit, timelimit);
	model->getEnv().set(GRB_IntParam_NumericFocus, 3);
	return 0;
}

int ilp_base::set_parameters()
{
	numericfocus = 3;
	//scaleflag = 0;
	gap = 0.0;
	nodes = 100000000;
	prepasses = 3;
	cuts = 2;
	heuristics = 0.0001;
	mipfocus = 2;
	//gomorypasses = 0;

	model->getEnv().set(GRB_IntParam_ScaleFlag, scaleflag);
	model->getEnv().set(GRB_DoubleParam_MIPGap, gap);
	model->getEnv().set(GRB_DoubleParam_NodeLimit, nodes);
	model->getEnv().set(GRB_IntParam_PrePasses, prepasses);
	model->getEnv().set(GRB_IntParam_Cuts, cuts);
	model->getEnv().set(GRB_IntParam_MIPFocus, mipfocus);
	model->getEnv().set(GRB_IntParam_GomoryPasses, gomorypasses);
	model->getEnv().set(GRB_DoubleParam_Heuristics, heuristics);
	return 0;
}

int ilp_base::check_connectivity()
{
	ugraph g(gr);
	int n = num_vertices(g);

	for(int i = 0; i < n; i++) add_edge(i, partners[i], g);

	vector<int> v(n);

	out_edge_iterator ei1, ei2;
	for(int i = 0; i < gene_pairs.size(); i++)
	{
		int g1 = gene_pairs[i].first;
		int g2 = gene_pairs[i].second;
		int g1h = gene_list[g1].first.first;
		int g1t = gene_list[g1].first.second;
		int g2h = gene_list[g2].first.first;
		int g2t = gene_list[g2].first.second;

		if(edge(g1h, g2h, gr).second == false) continue;
		if(edge(g1t, g2t, gr).second == false) continue;

		vector<int> e1h;
		vector<int> e1t;
		vector<int> e2h;
		vector<int> e2t;
		for(tie(ei1, ei2) = out_edges(g1h, gr); ei1 != ei2; ei1++) e1h.push_back(target(*ei1, gr));
		for(tie(ei1, ei2) = out_edges(g1t, gr); ei1 != ei2; ei1++) e1t.push_back(target(*ei1, gr));
		for(tie(ei1, ei2) = out_edges(g2h, gr); ei1 != ei2; ei1++) e2h.push_back(target(*ei1, gr));
		for(tie(ei1, ei2) = out_edges(g2t, gr); ei1 != ei2; ei1++) e2t.push_back(target(*ei1, gr));

		for(int j = i + 1; j < gene_pairs.size(); j++)
		{
			int g3 = gene_pairs[j].first;
			int g4 = gene_pairs[j].second;
			if(g3 == g1) continue;
			if(g4 == g2) continue;

			int g3h = gene_list[g3].first.first;
			int g3t = gene_list[g3].first.second;
			int g4h = gene_list[g4].first.first;
			int g4t = gene_list[g4].first.second;

			if(edge(g3h, g4h, gr).second == false) continue;
			if(edge(g3t, g4t, gr).second == false) continue;

			vector<int> e3h;
			vector<int> e3t;
			vector<int> e4h;
			vector<int> e4t;
			for(tie(ei1, ei2) = out_edges(g3h, gr); ei1 != ei2; ei1++) e3h.push_back(target(*ei1, gr));
			for(tie(ei1, ei2) = out_edges(g3t, gr); ei1 != ei2; ei1++) e3t.push_back(target(*ei1, gr));
			for(tie(ei1, ei2) = out_edges(g4h, gr); ei1 != ei2; ei1++) e4h.push_back(target(*ei1, gr));
			for(tie(ei1, ei2) = out_edges(g4t, gr); ei1 != ei2; ei1++) e4t.push_back(target(*ei1, gr));

			clear_vertex(g1h, g);
			clear_vertex(g1t, g);
			clear_vertex(g2h, g);
			clear_vertex(g2t, g);
			clear_vertex(g3h, g);
			clear_vertex(g3t, g);
			clear_vertex(g4h, g);
			clear_vertex(g4t, g);

			add_edge(g1h, g2h, g);
			add_edge(g1t, g2t, g);
			add_edge(g3h, g4h, g);
			add_edge(g3t, g4t, g);

			add_edge(g1h, partners[g1h], g);
			add_edge(g1t, partners[g1t], g);
			add_edge(g2h, partners[g2h], g);
			add_edge(g2t, partners[g2t], g);
			add_edge(g3h, partners[g3h], g);
			add_edge(g3t, partners[g3t], g);
			add_edge(g4h, partners[g4h], g);
			add_edge(g4t, partners[g4t], g);

			int c = connected_components(g, &v[0]);
			if(c == 1) continue;
			vector<int> cc;
			cc.assign(c, 0);
			printf("there are %6d components after fixing %5d and %5d pairs:", c, i, j);
			for(int k = 0; k < n; k++) cc[v[k]]++;
			for(int k = 0; k < c; k++) printf("%5d, ", cc[k]);
			printf("\n");

			for(int k = 0; k < e1h.size(); k++) add_edge(g1h, e1h[k], g);
			for(int k = 0; k < e1t.size(); k++) add_edge(g1t, e1t[k], g);
			for(int k = 0; k < e2h.size(); k++) add_edge(g2h, e2h[k], g);
			for(int k = 0; k < e2t.size(); k++) add_edge(g2t, e2t[k], g);
			for(int k = 0; k < e3h.size(); k++) add_edge(g3h, e3h[k], g);
			for(int k = 0; k < e3t.size(); k++) add_edge(g3t, e3t[k], g);
			for(int k = 0; k < e4h.size(); k++) add_edge(g4h, e4h[k], g);
			for(int k = 0; k < e4t.size(); k++) add_edge(g4t, e4t[k], g);
		}
	}

	return 0;
}

int ilp_base::print()
{
	// print graph
	printf("================= graph =================\n");
	for(int i = 0; i < num_vertices(gr); i++)
	{
		int g = gn_indices[i];
		char c;
		if(g == -1) c = '-';
		else if(gene_list[g].first.first == i) c = 'H';
		else if(gene_list[g].first.second == i) c = 'T';
		else assert(false);

		int gg = -1;
		if(g != -1) gg = gene_list[g].second;

		printf("vertex %4d [partner = %4d] [complement = %4d] [%4d:%4d:%c]: ", i, partners[i], complements[i], g, gg, c); 

		out_edge_iterator ei1, ei2;
		for(tie(ei1, ei2) = out_edges(i, gr); ei1 != ei2; ei1++)
		{
			int s = source(*ei1, gr);
			int t = target(*ei1, gr);
			assert(s == i);
			printf("%4d, ", t);
		}
		printf("\n");
	}

	printf("================= gene list =================\n");
	for(int i = 0; i < gene_list.size(); i++)
	{
		printf("gene: %4d, (%4d,%4d): %4d\n", i, gene_list[i].first.first, gene_list[i].first.second, gene_list[i].second);
	}

	return 0;
}
