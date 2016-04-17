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
	//genome(const genome_base & g, config * _conf);

	genome & operator= (const genome & g);

	virtual ~genome();

public:
	config * conf;

public:
	int find_max_gene();

	vector<int> build_gene_copy();
	int build_gene_map(vector< vector<gene*> > & list, gene* head, gene* tail);
	int build_gene_map(vector< vector<gene*> > & list);
	map<gene*, int> build_gene_indices();
	int build_gene_indices(map<gene*, int> & gi, map<int, gene*> & ig);

	int build_adjacencies(vector<adjacency> & s);

	int random_duplicate();

	int enumerate_operations(vector<operation*> & operations) const;

	int count_inversion();
	int count_duplication();
	int count_loss();

	int set_names(const string & prefix);
	int load(const string & file);
	int write(const string & file);

	int remove_tandem();
	bool remove_tandem_linear(chrm * ch, int p);

	set<gene*> merge_tandem();
	int merge_tandem(chrm * ch, set<gene*> & s);

	int statistic();
};

#endif
