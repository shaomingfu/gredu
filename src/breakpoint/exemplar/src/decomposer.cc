#include <fstream>
#include <iostream>
#include "decomposer.h"

decomposer::decomposer(genome * g1, genome * g2)
	: pbase(g1, g2)
{}

decomposer::~decomposer()
{}

int decomposer::solve()
{
	gm1->print();
	gm2->print();

	build();
	print_neato("dp.dot");
	return 0;
}

int decomposer::build()
{
	f2g.clear();
	g2f.clear();
	for(int i = 0;  i < ex1.size(); i++)
	{
		if(ex1[i] != NULL) continue;
		assert(ex2[i] == NULL);
		if(gf1[i].size() == 0) continue;
		if(gf2[i].size() == 0) continue;
		int n = f2g.size();
		f2g.insert(PI(i, n));
		g2f.insert(PI(n, i));
	}

	for(int k = 0; k < f2g.size(); k++) add_vertex(gr);

	for(int i = 0; i < vsa.size(); i++)
	{
		shadj sa = vsa[i];
		int f1 = (int)fabs(sa.x1->x);
		int f2 = (int)fabs(sa.x2->x);
		int g1 = f2g.find(f1) == f2g.end() ? -1 : f2g[f1];
		int g2 = f2g.find(f2) == f2g.end() ? -1 : f2g[f2];

		if(g1 != -1 && g2 != -1)
		{
			if(g1 != g2) add_edge(g1, g2, gr);
		}
		else if(g1 != -1 || g2 != -1)
		{
			int g = (g1 != -1) ? g1 : g2;
			vector<int> v1 = collect(sa.x1, sa.x2);
			vector<int> v2 = collect(sa.y1, sa.y2);
			for(int k = 0; k < v1.size(); k++)
			{
				int f = v1[k];
				if(f2g.find(f) == f2g.end()) continue;
				if(f2g[f] != g) add_edge(f2g[f], g, gr);
			}
			for(int k = 0; k < v2.size(); k++)
			{
				int f = v2[k];
				if(f2g.find(f) == f2g.end()) continue;
				if(f2g[f] != g) add_edge(f2g[f], g, gr);
			}
		}
		else
		{
			assert(g1 == -1 && g2 == -1);
			vector<int> v1 = collect(sa.x1, sa.x2);
			vector<int> v2 = collect(sa.y1, sa.y2);
			for(int k = 0; k < v1.size(); k++)
			{
				if(f2g.find(v1[k]) == f2g.end()) continue;
				for(int j = k + 1; j < v1.size(); j++)
				{
					if(f2g.find(v1[j]) == f2g.end()) continue;
					if(f2g[v1[j]] != f2g[v1[k]]) add_edge(f2g[v1[j]], f2g[v1[k]], gr);
				}
			}
			for(int k = 0; k < v2.size(); k++)
			{
				if(f2g.find(v2[k]) == f2g.end()) continue;
				for(int j = k + 1; j < v2.size(); j++)
				{
					if(f2g.find(v2[j]) == f2g.end()) continue;
					if(f2g[v2[j]] != f2g[v2[k]]) add_edge(f2g[v2[j]], f2g[v2[k]], gr);
				}
			}
		}
	}

	return 0;
}

vector<int> decomposer::collect(gene * x, gene * y)
{
	vector<int> v;
	PG p = sort(x, y);
	gene * g = p.first->b;
	while(g != p.second)
	{
		int f = (int)fabs(g->x);
		v.push_back(f);
		g = g->b;
	}
	return v;
}

int decomposer::print_neato(const string & file)
{
	ofstream fout(file.c_str());
	if(fout.fail()) return -1;
	edge_iterator ei1, ei2;
	fout<<"graph G {"<<endl;
	fout<<" graph [ overlap = false ];"<<endl;
	for(int i = 0; i < g2f.size(); i++)
	{
		int f = g2f[i];
		fout<<"v"<<i<<" [ label = "<<f<<", width = 0.1, height = 0.1 ];"<<endl;
	}

	for(tie(ei1, ei2) = edges(gr); ei1 != ei2; ei1++)
	{
		int s = source(*ei1, gr);
		int t = target(*ei1, gr);
		fout<<"v"<<s<<" -- "<<"v"<<t<<";"<<endl;
	}
	fout<<"}"<<endl;
	fout.close();
	return 0;
}
