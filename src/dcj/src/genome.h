#ifndef __GENOME_H__
#define __GENOME_H__

#include "genome_base.h"
#include "chrm.h"
#include "operation.h"
#include "config.h"
#include "adjacency.h"

#include <vector>
#include <map>
#include <set>

using namespace std;

class genome : public genome_base
{
public:
	genome(config * _conf);
	genome(const genome & g);

	genome & operator= (const genome & g);

	virtual ~genome();

public:
	config * conf;

public:
	int find_max_gene();

	vector<int> build_gene_copy();
	int build_gene_map(vector< vector<gene*> > & list);
	int build_gene_map(vector< vector<gene*> > & list, gene * head, gene * tail);

	int build_adjacencies(vector<adjacency> & s);

	int load(const string & file);
};


#endif
