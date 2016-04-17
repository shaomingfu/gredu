#ifndef __ASOLVER_H__
#define __ASOLVER_H__

#include "pbase.h"
#include "shseg.h"

class asolver : public pbase
{
public:
	asolver(config * _conf, genome * _gm1, genome * _gm2);
	virtual ~asolver();

public:
	ugraph ar;					// adjacency graph
	vector<int> adj;			// the adjacent extremity int the same adjacency
	vector<int> cpl;			// the complement extremity in the same gene

public:
	// build adjacency graph
	int build_adjacencies();
	int build_complements();
	int build_adjacency_graph();

	// utilities
	PG make_adjacency(int a, int b);
	PI make_adjacency(gene * a, gene * b);
	PI make_extremities(gene * g);
	bool adjacent(gene * x, gene * y);
	bool dcj_direction(int a, int b);
	bool check_availability(const PG & p);
	
	vector<PI> find_possible_adjacencies(int x, int y);
	bool shared(int x, int y);
	bool dominate(int x, int y);

	// resolve adjacencies
	int resolve_adjacencies();
	bool resolve_adjacency(int x, int y);
	int resolve_removed_genes();
	int resolve_removed_gene(int i);

	// identify shared adjacencies
	int identify_shared_adjacencies();
	bool identify_shared_adjacency(int x, int y);

	// identify duplicated circulars
	int identify_duplicated_circulars(const vector<PG> & dup);
	bool check_duplicated_circular(const PG & p);
	int remove_segment(const PG & p);

	// identify L-type
	int identify_ltypes();
	bool identify_ltype(int x, int y);

	// identify V-type
	int identify_vtypes();
	bool identify_vtype(int x, int y);
	
	// identify 8-type
	int identify_8types();
	bool identify_8type(int x, int y);

	// fix
	bool fix_extremity_pair(int a, int b);

	// print
	int print_genome(genome * gm);
	int print_adjacency_graph();
};

#endif
