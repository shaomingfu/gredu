#ifndef __GRED_H__
#define __GRED_H__

#include "distoptimizer.h"
#include "config.h"
#include "mygraph.h"

typedef map<gene*, gene*> MPG;

class gredo
{
public:
	gredo(double time);
	virtual ~gredo();
	
public:
	config * conf;
	genome * gm1;
	genome * gm2;

	ugraph gr;

public:
	int load_genomes(const string & file1, const string & file2);
	int solve(const string & file1, const string & file2);
	int write_mapping(const MPG & x2y, const string & file);
	int build_graph(const MPG & m);
};

#endif
