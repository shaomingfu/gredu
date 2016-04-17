#include "gredo.h"
#include "distoptimizer.h"

#include <fstream>
#include <algorithm>
#include <cstdio>

gredo::gredo(double time)
{
	conf = new config(time);
	gm1 = new genome(conf);
	gm2 = new genome(conf);
}

gredo::~gredo()
{
	if(conf != NULL) delete conf;
	if(gm1 != NULL) delete gm1;
	if(gm2 != NULL) delete gm2;
}

int gredo::load_genomes(const string & file1, const string & file2)
{
	gm1->load(file1);
	gm2->load(file2);
	int m1 = gm1->find_max_gene();
	int m2 = gm2->find_max_gene();
	conf->alphabet_size = m1 > m2 ? m1 : m2;
	gm1->alphabet_size = conf->alphabet_size;
	gm2->alphabet_size = conf->alphabet_size;
	return 0;
}

int gredo::solve(const string & file1, const string & file2)
{
	printf("load genomes ...\n");
	load_genomes(file1, file2);

	/*
	printf("genome1:\n");
	gm1->print();
	printf("genome2:\n");
	gm2->print();
	printf("alphabet size = %d\n", conf->alphabet_size);
	*/

	printf("run distoptimizer ...\n");
	distoptimizer d(conf, gm1, gm2);
	d.solve();

	write_mapping(d.x2y, "mapping");
	build_graph(d.x2y);
	return 0;
}

int gredo::write_mapping(const MPG & x2y, const string & file)
{
	ofstream fout(file.c_str());
	MPG::const_iterator it;
	for(it = x2y.begin(); it != x2y.end(); it++)
	{
		assert(it->first != NULL);
		assert(it->second!= NULL);
		if(it->first->x == 0) continue;
		fout<<it->first->s.c_str()<<" "<<it->second->s.c_str()<<endl;
	}
	fout.close();
	return 0;
}

int gredo::build_graph(const MPG & m)
{
	vector<adjacency> sx;
	vector<adjacency> sy;

	gm1->build_adjacencies(sx);
	gm2->build_adjacencies(sy);

	vector<int> copy1 = gm1->build_gene_copy();
	vector<int> copy2 = gm2->build_gene_copy();
	assert(copy1.size() == copy2.size());
	for(int i = 0; i < copy1.size(); i++) assert(copy1.at(i) == copy2.at(i));

	gr.clear();
	int n1 = sx.size() * 2;
	int n2 = sy.size() * 2;

	for(int i = 0; i < n1; i++) add_vertex(gr);
	for(int i = 0; i < n2; i++) add_vertex(gr);

	extremity ex, ey;
	for(int i = 0; i < sx.size(); i++)
	{
		ex = sx[i].e1;
		if(ex.g->x == 0) continue;
		assert(m.find(ex.g) != m.end());
		gene * g = m.find(ex.g)->second;
		ey = extremity(g, ex.b);

		int index = -1;
		for(int k = 0; k < sy.size(); k++)
		{
			if(sy[k].e1 == ey) index = n1 + k * 2;
			if(sy[k].e2 == ey) index = n1 + k * 2 + 1;
		}
		assert(index != -1);

		add_edge(i * 2, index, gr);
	}

	for(int i = 0; i < sx.size(); i++)
	{
		ex = sx[i].e2;
		if(ex.g->x == 0) continue;

		assert(m.find(ex.g) != m.end());
		gene * g = m.find(ex.g)->second;
		ey = extremity(g, ex.b);

		int index = -1;
		for(int k = 0; k < sy.size(); k++)
		{
			if(sy[k].e1 == ey) index = n1 + k * 2;
			if(sy[k].e2 == ey) index = n1 + k * 2 + 1;
		}
		assert(index != -1);

		add_edge(i * 2 + 1, index, gr);
	}

	assert(num_vertices(gr) % 2 == 0);
	for(int i = 0; i < num_vertices(gr) / 2; i++)
	{
		add_edge(i * 2, i * 2 + 1, gr);
	}

	return 0;
}
