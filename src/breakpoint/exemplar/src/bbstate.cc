#include "bbstate.h"
#include "adjacency.h"
#include <cmath>
#include <cassert>
#include <cstdio>

bbstate::bbstate(const vector<gene*> & _ig)
	: ig(_ig), distance(0)
{}

int bbstate::affect_candidate(int sa, int cd)
{
	if(cdm.find(sa) == cdm.end())
	{
		set<int> s;
		s.insert(cd);
		cdm.insert(pair<int, set<int> > (sa, s));
	}
	else
	{
		cdm[sa].insert(cd);
	}
	return 0;
}

int bbstate::add_pair(const PI & p)
{
	int inc = inc_distance(p);
	distance += inc;

	int x = p.first;
	int y = p.second;

	sx.insert(x);
	sy.insert(y);
	
	int f = (int)fabs(ig[x]->x);
	assert(f == (int)fabs(ig[y]->x));
	
	fx[f] = x;
	fy[f] = y;

	return 0;
}

int bbstate::add_pair(const PI & p, int inc)
{
	int x = p.first;
	int y = p.second;

	sx.insert(x);
	sy.insert(y);
	
	int f = (int)fabs(ig[x]->x);
	assert(f == (int)fabs(ig[y]->x));
	
	fx[f] = x;
	fy[f] = y;

	distance += inc;
	return 0;
}

int bbstate::remove_pair(const PI & p, int dec)
{
	int x = p.first;
	int y = p.second;

	assert(sx.find(x) != sx.end());
	assert(sy.find(y) != sy.end());
	sx.erase(x);
	sy.erase(y);
	
	int f = (int)fabs(ig[x]->x);
	assert(f == (int)fabs(ig[y]->x));
	
	assert(fx[f] != -1);
	assert(fy[f] != -1);

	fx[f] = -1;
	fy[f] = -1;

	distance -= dec;
	return 0;
}

int bbstate::remove_pair(const PI & p)
{
	int x = p.first;
	int y = p.second;

	assert(sx.find(x) != sx.end());
	assert(sy.find(y) != sy.end());
	sx.erase(x);
	sy.erase(y);
	
	int f = (int)fabs(ig[x]->x);
	assert(f == (int)fabs(ig[y]->x));
	
	assert(fx[f] != -1);
	assert(fy[f] != -1);

	fx[f] = -1;
	fy[f] = -1;

	int dec = inc_distance(p);

	distance -= dec;

	return 0;
}

int bbstate::inc_distance(const PI & p) const
{
	if(sx.size() == 0) return 0;

	int x = p.first;
	int y = p.second;

	gene * gx = ig[x];
	gene * gy = ig[y];

	set<int>::iterator ixu = sx.upper_bound(x);
	set<int>::iterator iyu = sy.upper_bound(y);

	int xu = (ixu == sx.end()) ? -1 : (*ixu);
	int yu = (iyu == sy.end()) ? -1 : (*iyu);
	int xl = (ixu == sx.begin()) ? -1 : (*(--ixu));
	int yl = (iyu == sy.begin()) ? -1 : (*(--iyu));

	int ans = 0;
	if(adjacent(x, xl, y, yl)) ans++;
	if(adjacent(x, xl, y, yu)) ans++;
	if(adjacent(x, xu, y, yl)) ans++;
	if(adjacent(x, xu, y, yu)) ans++;

	bool bx = false;
	if(xl != -1 && xu != -1)
	{
		int fl = (int)fabs(ig[xl]->x);
		int fu = (int)fabs(ig[xu]->x);
		int pl = fy[fl];
		int pu = fy[fu];
		if(adjacent(xl, xu, pl, pu)) bx = true;
	}

	bool by = false;
	bool bz = false;
	if(yl != -1 && yu != -1)
	{
		int fl = (int)fabs(ig[yl]->x);
		int fu = (int)fabs(ig[yu]->x);
		int pl = fx[fl];
		int pu = fx[fu];
		if(adjacent(pl, pu, yl, yu)) by = true;
		if(by && bx && pl == xl && pu == xu) bz = true;
		if(by && bx && pl == xu && pu == xl) bz = true;
	}

	if(bx == true) ans--;
	if(by == true) ans--;
	if(bz == true) ans++;
	
	return 1 - ans;
}

int bbstate::print() const
{
	int c = 0;
	for(int i = 0; i < cdb.size(); i++)
	{
		if(cdb[i] == true) c++;
	}

	int s = 0;
	for(int i = 0; i < sab.size(); i++)
	{
		if(sab[i] == true) s++;
	}


	printf("%d / %d families, %d breakpoints, %d / %d usable shared adjacency, %d / %d usable candidates\n", 
			(int)sx.size(), (int)(fx.size() - 1), distance, s, (int)sab.size(), c, (int)cdb.size());

	return 0;

	printf("gx: ");
	for(set<int>::iterator it = sx.begin(); it != sx.end(); it++)
	{
		int x = *it;
		int f = (int)fabs(ig[x]->x);
		printf("%3d(%3d) ", x, ig[x]->x);
	}
	printf("\n");

	printf("gy: ");
	for(set<int>::iterator it = sy.begin(); it != sy.end(); it++)
	{
		int y = *it;
		int f = (int)fabs(ig[y]->x);
		printf("%3d(%3d) ", y, ig[y]->x);
	}
	printf("\n");
	
	return 0;
}

bool bbstate::adjacent(int xl, int xu, int yl, int yu) const
{
	if(xl < 0 || xu < 0 || yl < 0 || yu < 0) return false;
	if(neighbor(sx, xl, xu) == false) return false;
	if(neighbor(sy, yl, yu) == false) return false;
	if(xl < xu && yl < yu && ig[xl]->x == ig[yl]->x && ig[xu]->x == ig[yu]->x) return true;
	if(xl > xu && yl > yu && ig[xl]->x == ig[yl]->x && ig[xu]->x == ig[yu]->x) return true;
	if(xl < xu && yl > yu && ig[xl]->x + ig[yl]->x == 0 && ig[xu]->x + ig[yu]->x == 0) return true;
	if(xl > xu && yl < yu && ig[xl]->x + ig[yl]->x == 0 && ig[xu]->x + ig[yu]->x == 0) return true;
	return false;
}

// return y in s s.t. y < x, return -1 if all elements are >= x
int find_less(const set<int> & s, int x)
{
	if(s.size() == 0) return -1;
	set<int>::const_iterator it = s.lower_bound(x);
	if(it == s.begin()) return -1;
	else return *(--it);
}

int find_less_equal(const set<int> & s, int x)
{
	if(s.size() == 0) return -1;
	if(s.find(x) != s.end()) return x;
	return find_less(s, x);
}

int find_greater(const set<int> & s, int x)
{
	if(s.size() == 0) return -1;
	set<int>::const_iterator it = s.upper_bound(x);
	if(it == s.end()) return -1;
	else return *it;
}

int find_greater_equal(const set<int> & s, int x)
{
	if(s.size() == 0) return -1;
	set<int>::const_iterator it = s.lower_bound(x);
	if(it == s.end()) return -1;
	else return *it;
}

bool neighbor(const set<int> & s, int x, int y)
{
	if(x == y) return false;
	int a = (x < y) ? x : y;
	int b = (x > y) ? x : y;
	if(find_greater(s, a) == b) return true;
	if(find_less(s, b) == a) return true;
	return false;
}
