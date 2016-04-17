#ifndef __PBASE0_H__
#define __PBASE0_H__

#include "pbase.h"
#include <vector>
#include <set>
#include <map>

using namespace std;

class pbase0 : public pbase
{
public:
	pbase0(config * conf, genome * g1, genome * g2, const MPG & _x2y, const MPG & _y2x);
	virtual ~pbase0();

public:
	vector<PG> pairs;
	map<PG, int> p2i;

	vector< set<int> > pct;

public:
	int build_gene_pairs();
	int add_gene_pair(const PG & p);
	int build_pairs_contact();
	int add_pair_contact(int g, int p);
};

#endif
