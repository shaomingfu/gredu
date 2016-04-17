#ifndef __PBASE_H__
#define __PBASE_H__

#include "genome.h"
#include "shadj.h"
#include "candidate.h"
#include <vector>
#include <set>
#include <map>

using namespace std;

class pbase
{
public:
	pbase(genome * g1, genome * g2);
	virtual ~pbase();

public:
	genome * gm1;
	genome * gm2;

	vector<gene*> ex1;
	vector<gene*> ex2;

	vector< vector<gene*> > gf1;
	vector< vector<gene*> > gf2;

	map<gene*, int> gi;
	vector<gene*> ig;

	vector<adjacency> va1;
	vector<adjacency> va2;
	vector<shadj> vsa;

	vector< set<int> > vsp;		// span list
	vector< set<int> > vct;		// contact list

	vector<candidate> cds;

public:
	int init();
	int build_gene_indices();
	int build_adjacencies(genome * gm, vector<adjacency> & va);
	int build_shared_adjacencies();
	int build_span_map();
	int build_candidates();
	bool build_candidate(const shadj & sa);
	bool jumpable(const candidate & cd);
	int fix_singletons();
	PG sort(gene * x, gene * y);
	int print_shared_adjacencies();
};

#endif
