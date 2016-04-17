#include "common.h"
#include <cmath>
#include <cfloat>
#include <cstdlib>
#include <cstdio>
#include <algorithm>

double log_add(double x, double y)
{
	if(x > y)
	{
		return x + log(1 + exp(y - x));
	}
	else
	{
		return y + log(1 + exp(x - y));
	}
}

int random(const vector<double> & ops)
{
	if(ops.size()==0) return -1;
	if(ops.size()==1) return 0;

	vector<double> vv;
	vv.assign(ops.size(), 0);

	double pre = 0;

	for(int i = 0; i < ops.size(); i++)
	{
		vv.at(i) = pre + ops.at(i);
		pre = vv.at(i);
	}

	double sum = vv.at(vv.size() - 1);

	for(int i = 0; i < vv.size(); i++)
	{
		//printf("vv.at(%d) = %6.3lf\n", i, vv.at(i));
		vv.at(i) = vv.at(i) / sum;
	}

	double r = (rand() % RAND_MAX) * 1.0 / RAND_MAX;

	for(int i = 0; i < vv.size(); i++)
	{
		if(r <= vv.at(i)) return i;
	}

	return vv.size() - 1;
}

int invert(vector<int> & v)
{
	std::reverse(v.begin(), v.end());
	
	for(int i = 0; i < v.size(); i++)
	{
		v.at(i) = 0 - v.at(i);
	}
	return 0;
}

