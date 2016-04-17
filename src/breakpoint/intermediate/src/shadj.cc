#include "shadj.h"
#include <cassert>
#include <algorithm>

shadj::shadj()
{}

shadj::shadj(gene * _x1, gene * _x2, gene * _y1, gene * _y2)
	: x1(_x1), x2(_x2), y1(_y1), y2(_y2)
{
	direction();
}

bool shadj::direction() const
{
	if(x1->x == y1->x && x2->x == y2->x) return true;
	else if(x1->x + y1->x == 0 && x2->x + y2->x == 0) return false;
	else assert(false);
}

bool shadj::adjacent() const
{
	if(x1->b != x2) return false;
	bool b = direction();
	if(b == true && y1->b == y2) return true;
	if(b == false && y2->b == y1) return true;
	return false;
}

bool shadj::conflict(const shadj & sa) const
{
	assert(adjacent());
	assert(sa.adjacent());
	if(x1 == sa.x1 && y1 == sa.y1) return true;
	if(x1 == sa.x2 && y1 == sa.y2) return true;
	if(x2 == sa.x1 && y2 == sa.y1) return true;
	if(x2 == sa.x2 && y2 == sa.y2) return true;
	if(x1 == sa.x1 || x1 == sa.x2) return false;
	if(x2 == sa.x1 || x2 == sa.x2) return false;
	if(y1 == sa.y1 || y1 == sa.y2) return false;
	if(y2 == sa.y1 || y2 == sa.y2) return false;
	return true;
}

bool shadj::coexist(const shadj & sa) const
{
	if(x1 == sa.x1 && y1 != sa.y1) return false;
	if(x1 == sa.x2 && y1 != sa.y2) return false;
	if(x2 == sa.x1 && y2 != sa.y1) return false;
	if(x2 == sa.x2 && y2 != sa.y2) return false;

	if(y1 == sa.y1 && x1 != sa.x1) return false;
	if(y1 == sa.y2 && x1 != sa.x2) return false;
	if(y2 == sa.y1 && x2 != sa.x1) return false;
	if(y2 == sa.y2 && x2 != sa.x2) return false;

	return true;
}

bool shadj::redundant() const
{
	vector<gene*> vx;
	vector<gene*> vy;
	vx = build_gene_list(PG(x1, x2));
	if(direction() == true) vy = build_gene_list(PG(y1, y2));
	else vy = build_gene_list(PG(y2, y1));

	set<int> sx;
	set<int> sy;
	for(int i = 1; i < vx.size() - 1; i++)
	{
		int x = vx[i]->x;
		if(sx.find(x) != sx.end()) continue;
		sx.insert(x);
	}

	for(int i = 1; i < vy.size() - 1; i++)
	{
		int x = vy[i]->x;
		if(direction() == false) x = 0 - x;
		if(sy.find(x) != sy.end()) continue;
		sy.insert(x);
	}

	set<int> s = set_intersection(sx, sy);
	if(s.size() >= 1) return true;
	else return false;
}

string shadj::print_string(map<gene*, int> & m) const
{
	char s[10240];
	sprintf(s, "(%6d[%6d], %6d[%6d]), (%6d[%6d], %6d[%6d])", m[x1], x1->x, m[x2], x2->x, m[y1], y1->x, m[y2], y2->x);
	return s;
}

vector<gene*> shadj::get_genes() const
{
	vector<gene*> v1 = build_gene_list(PG(x1, x2));
	vector<gene*> v2;
	if(direction() == true) v2 = build_gene_list(PG(y1, y2));
	else v2 = build_gene_list(PG(y2, y1));
	//printf("v1 = %5d, v2 = %5d\n", (int)v1.size(), (int)v2.size());
	v1.insert(v1.end(), v2.begin(), v2.end());
	return v1;
}

