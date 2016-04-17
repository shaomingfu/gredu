#include "candidate.h"
#include <algorithm>
#include <cassert>
#include <cstdio>
#include <cmath>

candidate::candidate()
{}

candidate::candidate(bool b, const vector<gene*> & x, const vector<gene*> & y)
	: xv(x), yv(y), d(b)
{
	assert(xv.size() == yv.size());
	
	/*
	for(int i = 0; i < x.size(); i++)
	{
		int f = (int)fabs(x[i]->x);
		assert(fs.find(f) == fs.end());
		fs.insert(f);
	}
	*/
}

bool candidate::operator<(const candidate & c) const
{
	if(xv.size() < c.xv.size()) return true;
	else return false;
}

/* 
 * comment because these are for examplars (handling gene familes)
bool candidate::conflict(const candidate & c) const
{
	vector<int> v;
	v.resize(size() + c.size());

	vector<int>::iterator it = set_intersection(fs.begin(), fs.end(), c.fs.begin(), c.fs.end(), v.begin());
	bool xy = (it - v.begin() >= 1) ? true : false;

	vector<int>::iterator xi = set_difference(fs.begin(), fs.end(), c.fs.begin(), c.fs.end(), v.begin());
	bool bx = (xi - v.begin() >= 1) ? true : false;

	vector<int>::iterator yi = set_difference(c.fs.begin(), c.fs.end(), fs.begin(), fs.end(), v.begin());
	bool by = (yi - v.begin() >= 1) ? true : false;

	if(bx && by && xy) return true;
	else return false;
}

int candidate::split(int p, candidate & c1, candidate & c2) const
{
	assert(p > 0 && p < size());
	vector<gene*> x1 (xv.begin(), xv.begin() + p);
	vector<gene*> y1 (yv.begin(), yv.begin() + p);
	c1 = candidate(d, x1, y1);

	vector<gene*> x2 (xv.begin() + p, xv.end());
	vector<gene*> y2 (yv.begin() + p, yv.end());
	c2 = candidate(d, x2, y2);

	return 0;
}
*/

int candidate::size() const
{
	return xv.size();
}

int candidate::print() const
{
	printf("[%d], ( ", d ? 1 : 0);
	for(int i = 0; i < xv.size(); i++) printf("%4d ", xv[i]->x);
	printf("), ( ");
	for(int i = 0; i < yv.size(); i++) printf("%4d ", yv[i]->x);
	printf(")\n");
	return 0;
}
