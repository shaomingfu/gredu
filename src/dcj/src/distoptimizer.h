#ifndef __DISTOPTIMIZER_H__
#define __DISTOPTIMIZER_H__

#include "mygraph.h"
#include "genome.h"
#include "config.h"

#include <vector>
#include <set>

using namespace std;

typedef pair<PI, int> PII;

class distoptimizer
{
public:
	distoptimizer(config * _conf, genome * _g1, genome * _g2);
	virtual ~distoptimizer();

public:
	config * conf;
	
	vector<adjacency> sx;
	vector<adjacency> sy;

	// ghost genes
	set<gene*> ghosts;

	// result mapping
	map<gene*, gene*> x2y;
	map<gene*, gene*> y2x;

	// input for ILP
	vector<PII> gene_list;
	vector<int> partners;

	int simplified_cycles;

	ugraph gr;				// simplified graph
	vector<int> b2a;		// index map from simplified to raw graph
	vector<int> a2b;		// index map from raw to simplified graph

	ugraph bg;				// final breakpoint graph

public:
	int solve();

	int check_components();

	int build_adjacency_sets(genome * gm1, genome * gm2);

	int check_edge(const extremity & ex, const extremity & ey);
	int build_graph();
	int simplify_graph();
	int build_gene_list();

	int build_mapping(const set<PI> & pairs);
	int add_gene_pair(gene * gx, gene * gy);
	int add_extremity_pair(int x, int y);
};

#endif
