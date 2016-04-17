#include "trier.h"
#include "simulator.h"
#include "exemplar.h"

#include <fstream>
#include <algorithm>
#include <cstdio>

trier::trier(const char * file)
{
	conf = new config(file);
	gm1 = new genome(conf);
	gm2 = new genome(conf);
}

trier::trier(double timelimit)
{
	conf = new config();
	conf->ilp_timelimit = timelimit;

	gm1 = new genome(conf);
	gm2 = new genome(conf);
}

trier::~trier()
{
	if(conf != NULL) delete conf;
	if(gm1 != NULL) delete gm1;
	if(gm2 != NULL) delete gm2;
}

int trier::load_genomes(const string & file1, const string & file2)
{
	gm1->load(file1);
	gm2->load(file2);
	int m1 = gm1->find_max_gene();
	int m2 = gm2->find_max_gene();
	conf->alphabet_size = m1 > m2 ? m1 : m2;
	return 0;
}

int trier::simulate_genomes()
{
	simulator sm(conf, gm1, gm2);
	sm.simulate();
	s2t = sm.s2t;
	t2s = sm.t2s;
	gm1->write("gm1", "X");
	gm2->write("gm2", "Y");

	int sa = gm1->sadist(*gm2, s2t);

	set<gene*> ss;
	set<gene*> tt;
	map<gene*, gene*>::iterator it;
	for(it = s2t.begin(); it != s2t.end(); it++)
	{
		ss.insert(it->first);
		tt.insert(it->second);
	}

	gm1->printp("XXX");
	gm2->printp("XXX");
	printf("\n");

	gm1->print();
	printf("---\n");
	gm2->print();
	printf("\n");

	printf("oracle genomes ...\n");
	gm1->prints(ss);
	printf("---\n");
	gm2->prints(tt);

	printf("oracle  shared adjacencies= %5d\n\n", sa);

	return 0;
}

int trier::solve()
{
	printf("simulate genomes ...\n");
	simulate_genomes();

	printf("run exemplar ...\n");
	exemplar d(conf, gm1, gm2);
	d.solve();

	x2y = d.x2y;
	y2x = d.y2x;

	/*
	printf("\n");
	gm1->print();
	printf("---\n");
	gm2->print();
	*/

	printf("optimal shared adjacencies = %5d, remaining = %5d\n", d.nshadj, gm1->sadist(*gm2, x2y));

	printf("\nevaluate ...\n");
	evaluate();
	//write_mapping(x2y, "mapping");

	return 0;
}

int trier::solve(const string & file1, const string & file2)
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

	printf("run exemplar ...\n");
	exemplar d(conf, gm1, gm2);
	d.solve();

	x2y = d.x2y;
	y2x = d.y2x;

	/*
	printf("\n");
	gm1->print();
	printf("---\n");
	gm2->print();
	*/

	printf("optimal shared adjacencies = %5d, remaining = %5d\n", d.nshadj, gm1->sadist(*gm2, x2y));

	//evaluate();
	write_mapping(x2y, "mapping");

	return 0;
}

int trier::evaluate()
{
	//assert(ex1.size() == s2t.size());
	//assert(ex2.size() == t2s.size());
	map<gene*, gene*>::iterator it;
	int s = 0;
	int s1 = 0;
	int s2 = 0;
	for(it = x2y.begin(); it != x2y.end(); it++)
	{
		if(s2t.find(it->first) != s2t.end() && s2t[it->first] == it->second) s++;
	}
	for(it = x2y.begin(); it != x2y.end(); it++)
	{
		if(s2t.find(it->first) != s2t.end()) s1++;
	}
	for(it = y2x.begin(); it != y2x.end(); it++)
	{
		if(t2s.find(it->first) != t2s.end()) s2++;
	}

	int n = x2y.size();
	printf("genomes: %5d out of %5d     pairs are correct. ratio = %8.3lf\n", s,  n, 1.0 * s  / n);
	printf("genome1: %5d out of %5d exemplars are correct. ratio = %8.3lf\n", s1, n, 1.0 * s1 / n);
	printf("genome2: %5d out of %5d exemplars are correct. ratio = %8.3lf\n", s2, n, 1.0 * s2 / n);

	return 0;
}

int trier::write_mapping(const map<gene*, gene*> & m, const string & file)
{
	ofstream fout(file.c_str());
	map<gene*, gene*>::const_iterator it;
	for(it = m.begin(); it != m.end(); it++)
	{
		assert(it->first != NULL);
		assert(it->second!= NULL);
		if(it->first->x == 0) continue;
		fout<<it->first->s.c_str()<<" "<<it->second->s.c_str()<<endl;
	}
	fout.close();
	return 0;
}
