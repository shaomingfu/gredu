#include "common.h"
#include <cassert>
#include <algorithm>

set<int> set_intersection(const set<int> & x, const set<int> & y)
{
	vector<int> z(x.size() + y.size());
	vector<int>::iterator end = std::set_intersection(x.begin(), x.end(), y.begin(), y.end(), z.begin());
	return set<int>(z.begin(), end);
}

set<int> set_union(const set<int> & x, const set<int> & y)
{
	vector<int> z(x.size() + y.size());
	vector<int>::iterator end = std::set_union(x.begin(), x.end(), y.begin(), y.end(), z.begin());
	return set<int>(z.begin(), end);
}
