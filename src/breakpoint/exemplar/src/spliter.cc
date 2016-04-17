#include "spliter.h"
#include <algorithm>
#include <cmath>

spliter::spliter(genome * gm1, genome * gm2)
	: pbase(gm1, gm2)
{
	split();
}

spliter::~spliter()
{}

int spliter::split()
{
	std::sort(cds.begin(), cds.end());
	//print();

	add();
	
	//print();

	for(int i = 0; i < cps.size(); i++)
	{
		for(int j = 0; j < cps.size(); j++)
		{
			bool b = cps[i].conflict(cps[j]);
			assert(b == false);
		}
	}

	return 0;
}

int spliter::add()
{
	cps.clear();
	if(cds.size() <= 0) return 0;

	cps.push_back(cds[0]);
	for(int k = 1; k < cds.size(); k++)
	{
		assert(cds[k].size() >= 2);
		add(cds[k]);
	}

	return 0;
}

int spliter::add(const candidate & y)
{
	for(int k = 0; k < cps.size(); k++)
	{
		if(cps[k].conflict(y) == false) continue;
		vector<candidate> v = compare(cps[k], y);
		for(int i = 0; i < v.size(); i++) add(v[i]);
		return 0;
	}

	if(y.size() >= 2) cps.push_back(y);
	return 0;
}

vector<candidate> spliter::compare(const candidate & x, const candidate & y) const
{
	vector<candidate> v;
	if(x.conflict(y) == false)
	{
		v.push_back(y);
		return v;
	}

	candidate c1, c2;
	for(int p = y.size() - 1; p >= 1; p--)
	{
		y.split(p, c1, c2);
		if(x.conflict(c1) == false) break;
	}

	v = compare(x, c2);
	v.push_back(c1);

	return v;
}

int spliter::print() const
{
	printf(" ------------ original candidates -----------\n");
	for(int i = 0; i < cds.size(); i++)
	{
		printf("%5d:", i + 1);
		cds[i].print();	
	}

	printf(" ------------ splitted candidates -----------\n");
	for(int i = 0; i < cps.size(); i++)
	{
		printf("%5d:", i + 1);
		cps[i].print();	
	}

	return 0;
}

