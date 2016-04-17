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

typedef map<gene*, gene*> MPG;

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
	map<gene*, int> build_gene_indices();
	int build_gene_map(vector< vector<gene*> > & list);
	set<gene*> build_gene_set();

	int random_duplicate();

	int build_adjacencies(vector<adjacency> & s) const;
	int build_adjacencies(vector<adjacency> & s, const set<gene*> & st) const;

	// enumerate all the possible operations
	int enumerate_operations(vector<operation*> & operations) const;

	int count_inversion();
	int count_duplication();
	int count_loss();

	int load(const string & file);
	int write(const string & file, const string & prefix);
	int size();

	int remove_tandem();
	bool remove_tandem_linear(chrm * ch, int p);
	bool remove_tandem_linear_triple(chrm * ch);

	int sadist(const genome & gm, const map<gene*, gene*> & m);

	int prints(const set<gene*> & s) const;
	int printp(const string & p) const;
	int printc(int cutoff) const;
	int printi(map<gene*, int> & m) const;
};


#endif
