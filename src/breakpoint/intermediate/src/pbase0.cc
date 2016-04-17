#include "pbase0.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <algorithm>

using namespace std;

pbase0::pbase0(config * conf, genome * g1, genome * g2, const MPG & _x2y, const MPG & _y2x)
	: pbase(conf, g1, g2, _x2y, _y2x)
{
	build_gene_pairs();
	build_pairs_contact();
}

pbase0::~pbase0()
{}

int pbase0::add_gene_pair(const PG & p)
{
	if(p2i.find(p) != p2i.end()) return p2i[p];
	int n = p2i.size();
	p2i.insert(pair<PG, int>(p, n));
	return n;
}

int pbase0::build_gene_pairs()
{
	for(MPG::iterator it = x2y.begin(); it != x2y.end(); it++)
	{
		gene * x = it->first;
		gene * y = it->second;
		add_gene_pair(PG(x, y));
	}

	for(int k = 0; k < vsa.size(); k++)
	{
		shadj & sa = vsa[k];
		PG p1(sa.x1, sa.y1);
		PG p2(sa.x2, sa.y2);
		add_gene_pair(p1);
		add_gene_pair(p2);
	}
	return 0;
}

int pbase0::add_pair_contact(int g, int p)
{
	assert(g >= 0 && g < pct.size());
	set<int> & s = pct[g];
	if(s.find(p) != s.end()) return 0;
	s.insert(p);
	return 0;
}

int pbase0::build_pairs_contact()
{
	pct.clear();
	pct.resize(gi.size());
	map<PG, int>::iterator it;
	for(it = p2i.begin(); it != p2i.end(); it++)
	{
		PG p = it->first;
		int k = it->second;
		int gx = gi[p.first];
		int gy = gi[p.second];
		add_pair_contact(gx, k);
		add_pair_contact(gy, k);
	}
	return 0;
}
