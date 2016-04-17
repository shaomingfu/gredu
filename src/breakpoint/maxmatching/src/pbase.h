#ifndef __PBASE_H__
#define __PBASE_H__

#include "config.h"
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
	pbase(config * c, genome * g1, genome * g2, const MPG & _x2y, const MPG & _y2x);
	virtual ~pbase();

public:
	config * conf;
	genome * gm1;
	genome * gm2;
	vector<int> nm;

	MPG x2y;
	MPG y2x;
	vector<int> vf;

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
	int run();

	int init_num_matching();
	int init_fixed_families();
	int build_gene_indices();
	int build_adjacencies(genome * gm, vector<adjacency> & va);
	int build_shared_adjacencies();
	int build_contact_map();
	int build_span_map();
	int build_candidates();
	bool build_candidate(const shadj & sa);
	PG sort(gene * x, gene * y);

	bool check_redundance(const shadj & sa, const vector<bool> & s);
	int remove_redundant_shared_adjacencies();
	int remove_redundant_genes();

	int print_families();
	int print_shared_adjacencies();
	int print_shared_adjacencies(const set<int> & s);
	//bool jumpable(const candidate & cd);
	
	bool innocent(const shadj & sa, const vector<gene*> & xv, const vector<gene*> & yv);

	bool is_fixed_gene(gene * x);
	bool is_fixed_pair(gene * x, gene * y);
	bool linkable(gene * x, gene * y);
	int add_fixed_pair(gene * x, gene * y);

	int fix_singletons();

	int locate_shared_adjacency(gene * x1, gene * x2, gene * y1, gene * y2);
	vector<int> locate_shared_adjacencies(const candidate & c);
	vector<int> locate_shared_adjacencies(const vector<gene*> & xv, const vector<gene*> & yv);

	set<int> build_pair_contact(gene * x, gene * y);
	set<int> build_span_intersection(const vector<gene*> & v);
	set<int> build_span_union(const vector<gene*> & v);
	set<int> build_contact_intersection(const vector<gene*> & v);
	set<int> build_contact_union(const vector<gene*> & v);
};

 #endif
