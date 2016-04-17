#ifndef __GFAMILY_H__
#define __GFAMILY_H__

#include <vector>

using namespace std;

typedef pair<int, int> PI;

class gpair
{
public:
	gpair(const PI & _p, int _s);
	gpair(const PI & _p);
	bool operator<(const gpair & g) const;

public:
	PI p;
	int s;
};

class gfamily
{
public:
	int f;
	vector<gpair> vp;

public:
	gfamily(int _f);
	int size() const;
	bool operator<(const gfamily & gf) const;
	int sort();
	int print();
};

#endif
