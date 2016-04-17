#include "shadj.h"
#include <cassert>

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

string shadj::print_string(map<gene*, int> & m) const
{
	char s[10240];
	sprintf(s, "(%7d[%7d], %7d[%7d]), (%7d[%7d], %7d[%7d])", m[x1], x1->x, m[x2], x2->x, m[y1], y1->x, m[y2], y2->x);
	return s;
}
