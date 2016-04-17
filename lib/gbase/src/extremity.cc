#include "extremity.h"
#include <cmath>
#include <cassert>
#include <cstdio>

using namespace std;

extremity::extremity()
{
	g = NULL;
	b = true;
}

extremity::extremity(gene * _g, bool _b)
{
	g = _g;
	b = _b;
}

extremity::extremity(const extremity & e)
{
	g = e.g;
	b = e.b;
}

bool extremity::operator==(const extremity & e) const
{
	if(g->x == 0 && g == e.g) return true;
	if(g == e.g && b == e.b) return true;
	return false;
}

bool extremity::operator!=(const extremity & e) const
{
	if((*this) == e) return false;
	return true;
}

bool extremity::forward() const
{
	if(g->x == 0)
	{
		if(g->b == NULL) return false;
		else if(g->a == NULL) return true;
		else assert(1 == 0);
	}
	if(g->x > 0 && b == true) return true;
	if(g->x < 0 && b == false) return true;
	return false;
}

bool extremity::complement(const extremity & e) const
{
	if(g->x == 0) return false;
	if(g == e.g && b != e.b) return true;
	return false;
}

int extremity::label() const
{
	if(forward()) return g->x;
	else return 0 - g->x;
}

string extremity::label_tex() const
{
	if(g->x == 0) return "00";
	else
	{
		char c;
		if(b == true) c = 'h';
		else c = 't';

		char buf[1024];
		sprintf(buf, "%d^%d_%c", (int)(fabs(g->x)), gene::get_gene_copy_index(g), c);
		return buf;
	}
}

string extremity::label_str() const
{
	char c;
	if(b == true) c = 'H';
	else c = 'T';

	char buf[1024];
	//sprintf(buf, "%d%c", (int)(fabs(g->x)), c);
	int d = 0;
	if(g->x != 0) d = gene::get_gene_copy_index(g);
	sprintf(buf, "%d%c%d", (int)(fabs(g->x)), c, d);
	return buf;
}

int extremity::weak_compare(const extremity & e) const
{
	if((int)fabs(e.g->x) != (int)fabs(g->x)) return 0;
	if(g->x == 0) return 1;
	if(e.b == b) return 1;
	return 0;
}

