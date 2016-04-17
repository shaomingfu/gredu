#include "gfamily.h"
#include <cassert>
#include <algorithm>
#include <cstdio>

gpair::gpair(const PI & _p, int _s)
	: p(_p), s(_s)
{}

gpair::gpair(const PI & _p)
	: p(_p)
{}

bool gpair::operator<(const gpair & g) const
{
	if(s < g.s) return true;
	else return false;
}

gfamily::gfamily(int _f)
	:f(_f)
{}

bool gfamily::operator<(const gfamily & gf) const
{
	if(vp.size() < gf.vp.size()) return true;
	return false;
}

int gfamily::size() const
{
	return vp.size();
}

int gfamily::sort()
{
	std::sort(vp.begin(), vp.end());
	return 0;
}

int gfamily::print()
{
	for(int i = 0; i < vp.size(); i++)
	{
		printf("pair (%4d, %4d), score = %4d\n", vp[i].p.first, vp[i].p.second, vp[i].s);
	}
	return 0;
}
