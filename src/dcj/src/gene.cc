#include "gene.h"
#include "chrm.h"

#include <cmath>
#include <cstddef>
#include <map>
#include <vector>

gene::gene(int _x)
{
	x = _x;
	a = NULL;
	b = NULL;
	ch = NULL;
}

gene::gene(int _x, chrm * _ch, const string & _s)
{
	x = _x;
	a = NULL;
	b = NULL;
	ch = _ch;
	s = _s;
}

gene::gene(int _x, chrm * _ch)
{
	x = _x;
	a = NULL;
	b = NULL;
	ch = _ch;
}

/*
 * dangerous!!!
gene::gene(const gene & g)
{
	x = g.x;
	a = g.a;
	b = g.b;
	ch = g.ch;
}
*/

int gene::reverse()
{
	x = 0 - x;
	gene * t = a;
	a = b;
	b = t;
	return 0;
}

int gene::get_gene_copy_index(gene * g)
{
	typedef vector<gene*> VG;
	typedef map<int, VG> MIVG;

	static MIVG m;
	int x = (int)(fabs(g->x));

	MIVG::iterator it = m.find(x);

	if(it == m.end())
	{
		vector<gene*> v;
		v.push_back(g);
		m.insert(pair<int, VG>(x, v));
		return 0;
	}
	else
	{
		for(int i = 0; i < it->second.size(); i++)
		{
			if(it->second.at(i) == g) return i;
		}
		it->second.push_back(g);
		return it->second.size() - 1;
	}
}

vector<gene*> build_gene_list(const PG& pg)
{
	gene * p = pg.first;
	vector<gene*> v;
	while(true)
	{
		v.push_back(p);	
		if(p == pg.second) break;
		p = p->b;
	}
	return v;
}

