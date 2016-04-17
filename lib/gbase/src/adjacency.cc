#include "adjacency.h"
#include "common.h"
#include "chrm.h"

#include <cassert>
#include <cstdio>

adjacency::adjacency(const adjacency & a)
{
	e1 = a.e1;
	e2 = a.e2;
}

adjacency::adjacency(const extremity & _e1, const extremity & _e2)
{
	e1 = _e1;
	e2 = _e2;
}

adjacency & adjacency::operator=(const adjacency & a)
{
	e1 = a.e1;
	e2 = a.e2;
	return (*this);
}

adjacency::adjacency()
{}

adjacency::adjacency(const PG & pg)
{
	e1.g = pg.first;
	e2.g = pg.second;

	if(e1.g->x >= 0) e1.b = true;
	else e1.b = false;

	if(e2.g->x >= 0) e2.b = false;
	else e2.b = true;
}

bool adjacency::operator!=(const adjacency & a) const
{
	if( (*this) == a ) return false;
	else return true;
}

bool adjacency::operator==(const adjacency & a) const
{
	if(e1 == a.e1 && e2 == a.e2) return true;
	if(e1 == a.e2 && e2 == a.e1) return true;
	return false;
}

bool adjacency::available() const
{
	if(e1.g->b == e2.g && e1.forward() == true && e2.forward() == false) return true;
	if(e2.g->b == e1.g && e2.forward() == true && e1.forward() == false) return true;
	return false;
}

int adjacency::print() const
{
	printf("(%3d,%3d)", e1.label(), 0 - e2.label());
	return 0;
}

int adjacency::exchange()
{
	extremity e = e1;
	e1 = e2;
	e2 = e;
	return 0;
}

string adjacency::label_tex() const
{
	char buf[1024];
	sprintf(buf, "%s%s", e1.label_tex().c_str(), e2.label_tex().c_str());
	//sprintf(buf, "(%d,%d)", e1.label(), 0 - e2.label());
	return string(buf);
}

string adjacency::label_str() const
{
	char buf[1024];
	sprintf(buf, "(%s,%s)", e1.label_str().c_str(), e2.label_str().c_str());
	//sprintf(buf, "(%d,%d)", e1.label(), 0 - e2.label());
	return string(buf);
}

bool adjacency::direction() const
{
	bool b1 = e1.forward();
	bool b2 = e2.forward();

	assert(b1 != b2);
	return b1;
}

int adjacency::strong_compare(const adjacency & a) const
{
	if(e1 == a.e1 && e2 == a.e2) return 2;
	if(e1 == a.e2 && e2 == a.e1) return 2;
	if(e1 == a.e1 || e1 == a.e2) return 1;
	if(e2 == a.e1 || e2 == a.e2) return 1;
	return 0;
}

int adjacency::weak_compare(const adjacency & a) const
{
	int x1 = e1.weak_compare(a.e1) + e2.weak_compare(a.e2);
	int x2 = e1.weak_compare(a.e2) + e2.weak_compare(a.e1);
	return x1 > x2 ? x1 : x2;
}

int adjacency::sort_adjacencies(vector<adjacency> & s)
{
	if(s.size() <= 1) return 0;

	vector<int> list;
	vector<bool> mask;
	mask.assign(s.size(), false);

	while(true)
	{
		if(list.size() == s.size()) break;

		// finding the starting point which has no predecessor and no successor
		int p = -1;
		for(int i = 0; i < s.size(); i++)
		{
			if(mask.at(i) == true) continue;
			p = i;

			bool b1 = false;
			bool b2 = false;
			for(int j = 0; j < s.size(); j++)
			{
				if(i == j) continue;
				if(mask.at(j) == true) continue; 
				if(s.at(i).e1.complement(s.at(j).e1)) b1 = true;
				if(s.at(i).e1.complement(s.at(j).e2)) b1 = true;
				if(s.at(i).e2.complement(s.at(j).e1)) b2 = true;
				if(s.at(i).e2.complement(s.at(j).e2)) b2 = true;
			}

			if(b1 == true && b2 == true) continue;

			if(b1 == true) s.at(p).exchange();
			break;
		}

		list.push_back(p);
		mask.at(p) = true;

		while(true)
		{
			bool found = false;
			for(int i = 0; i < s.size(); i++)
			{
				if(mask.at(i) == true) continue;

				if(s.at(p).e2.complement(s.at(i).e1))
				{
					p = i;
					found = true;
					break;
				}

				if(s.at(p).e2.complement(s.at(i).e2))
				{
					p = i;
					s.at(p).exchange();
					found = true;
					break;
				}
			}

			if(found == false) break;
			list.push_back(p);
			assert(mask.at(p) == false);
			mask.at(p) = true;
		}
	}

	vector<adjacency> ss = s;
	for(int i = 0; i < s.size(); i++)
	{
		s.at(i) = ss.at(list.at(i));
	}

	return 0;
}

adjacency adjacency::diff(const adjacency & a, const adjacency & o) const
{
	extremity x1, x2;
	if(e1 == o.e1 || e1 == o.e2) x1 = e2;
	else if(e2 == o.e1 || e2 == o.e2) x1 = e1;
	else assert(1 == 0);

	if(a.e1 == o.e1 || a.e1 == o.e2) x2 = a.e2;
	else if(a.e2 == o.e1 || a.e2 == o.e2) x2 = a.e1;
	else assert(1 == 0);

	return adjacency(x1, x2);
}

extremity adjacency::diff(const adjacency & o) const
{
	assert(strong_compare(o) == 1);
	if(e1 == o.e1 || e1 == o.e2) return e2;
	else if(e2 == o.e1 || e2 == o.e2) return e1;
	else assert(1 == 0);
}

extremity adjacency::same(const adjacency & o) const
{
	assert(strong_compare(o) == 1);
	if(e1 == o.e1 || e1 == o.e2) return e1;
	else if(e2 == o.e1 || e2 == o.e2) return e2;
	else assert(1 == 0);
}

vector<adjacency> build_adjacency_list(const PG & pg)
{
	gene * p = pg.first;
	vector<adjacency> v;
	while(true)
	{
		if(p == pg.second) break;
		v.push_back(adjacency(PG(p,p->b)));
		p = p->b;
	}
	return v;
}
