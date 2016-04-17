#include "trier.h"
#include "simulator.h"
#include "psolver.h"
#include "lpsolver.h"
#include "lpsolver2.h"

#include <fstream>
#include <algorithm>
#include <cstdio>

trier::trier(const char * file)
{
	conf = new config(file);
	gm1 = new genome(conf);
	gm2 = new genome(conf);
}

trier::trier(double limit)
{
	conf = new config(limit);
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
	gm1->set_names("X");
	gm2->set_names("Y");
	gm1->write("gm1");
	gm2->write("gm2");

	write_mapping(sm.s2t, "map");

	return 0;
}

int trier::solve()
{
	printf("simulate genomes ...\n");
	simulate_genomes();
	printf("\n");
	return 0;
}

int trier::solve(const string & file1, const string & file2)
{
	printf("load genomes ...\n");
	load_genomes(file1, file2);
	printf("alphabet size = %d\n", conf->alphabet_size);

	printf("run lpsolver ...\n");
	lpsolver lp(conf, gm1, gm2);
	lp.solve();
	write_mapping(lp.x2y, "mapping");
	//write_mapping(lp.px2y, "fixed");

	return 0;
}

int trier::write_mapping(const MPG & x2y, const string & file)
{
	ofstream fout(file.c_str());
	MPG::const_iterator it;
	for(it = x2y.begin(); it != x2y.end(); it++)
	{
		assert(it->first != NULL);
		//assert(it->second!= NULL);
		if(it->first->x == 0) continue;
		if(it->second == NULL) continue;
		fout<<it->first->s.c_str()<<" "<<it->second->s.c_str()<<endl;
	}
	fout.close();
	return 0;
}
