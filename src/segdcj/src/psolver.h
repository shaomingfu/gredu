#ifndef __PSOLVER_H__
#define __PSOLVER_H__

#include "asolver.h"
#include "bsolver.h"

typedef vector< vector<int> > VVI;
typedef vector< vector<gene*> > VVG;

class psolver: public bsolver
{
public:
	psolver(config * _conf, genome * _gm1, genome * _gm2);
	virtual ~psolver();

public:
	// prepared for lpsolver2
	VVI gc1;			// gene families1
	VVI gc2;			// gene families2
	VVI vdup;			// duplicons
	vector<PI> ge;		// gene extremities
	vector<int> gf;		// gene families
	vector<PI> vip;		// inference pairs

public:
	int build_gene_list();
	int build_gene_families();
	int transform_duplicons();

	// inference pairs
	int build_inference_pairs();
	bool add_inference_pair(int s, int t, gene * x, gene * y);

	map<gene*, gene*> transform_mapping(const map<int, int> & m);

	int check_families();
};

#endif
