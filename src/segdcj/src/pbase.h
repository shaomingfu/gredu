#ifndef __PBASE_H__
#define __PBASE_H__

#include "mygraph.h"
#include "genome.h"
#include "config.h"

typedef vector< vector<gene*> > VVG;

class pbase
{
public:
	pbase(config * _conf, genome * _gm1, genome * _gm2);
	virtual ~pbase();

public:
	genome * gm1;
	genome * gm2;

	config * conf;

	map<gene*, int> gi1;				// all genes, including 0s
	map<gene*, int> gi2;				// all genes, including 0s
	map<int, gene*> ig1;				// gene map 
	map<int, gene*> ig2;				// gene map 

	ugraph pr;							// gene graph
	map<edge_descriptor, int> ei;		// edge map
	map<int, edge_descriptor> ie;		// edge map

	VVG gf1;							// gene families1
	VVG gf2;							// gene families2
	vector<PG> dup1;					// duplicons in genome1
	vector<PG> dup2;					// duplicons in genome2

	map<gene*, gene*> px2y;
	map<gene*, gene*> py2x;


public:
	// make the 0s equal in the two genomes
	int add_ghost_chrms();

	// remove tandem arrayed genes
	int remove_tandem();

	// gene graph
	int build_gene_indices();
	int build_edge_map();
	int build_gene_graph();
	int build_gene_graph(const vector<PG> & v);
	int store_gene_graph(vector<PG> & v);
	
	// utilities
	int get_index(gene * g);
	gene * get_gene(int i);
	int get_equal_size(int g);
	int get_equal_size(gene * g);
	int get_copy_number(int g);
	int get_copy_number(gene * g);
	bool check_all_duplicated(const PG & p);

	// build duplicons
	int build_duplicons(const VVG & gf, vector<PG> & dup);
	bool check_duplicon(const PG & p);

	// build prefixed mappings
	int build_fixed_gene_pairs();

	// fixing gene
	bool is_gene_fixed(int x);
	bool is_gene_pair_fixed(gene * x, gene * y);
	bool fix_gene_pair(gene * x, gene * y);
	bool fix_gene_pair(int x, int y);

	// statistics
	int statistic_degrees();
	int statistic_family(const vector<gene*> & v, map<gene*, int> & gi);
	int statistic_forks();
	int statistic_shared_adjacencies();
	int statistic_fixed_gene_pairs();

	// print
	int print_gene_graph();
	string print_triple(gene * x, map<gene*, int> & gi, int n);
};

#endif
