#ifndef __SHADJ__
#define __SHADJ__

#include "common.h"
#include "gene.h"
#include <set>
#include <map>
#include <cstdio>
#include <string>

using namespace std;

typedef map<gene*, gene*> MPG;

class shadj
{
public:
	shadj();
	shadj(gene * _x1, gene * _x2, gene * _y1, gene * _y2);

public:
	gene * x1;
	gene * x2;
	gene * y1;
	gene * y2;

public:
	bool direction() const;
	bool adjacent() const;
	bool conflict(const shadj & sa) const;
	bool coexist(const shadj & sa) const;
	bool redundant() const;
	vector<gene*> get_genes() const;
	string print_string(map<gene*, int> & m) const;
};

#endif
