#ifndef __TRIER_H__
#define __TRIER_H__

#include "intermediate.h"
#include "config.h"

class trier
{
public:
	trier(double timelimit);
	trier(const char * file);
	virtual ~trier();
	
public:
	config * conf;
	genome * gm1;
	genome * gm2;

	map<gene*, gene*> x2y;
	map<gene*, gene*> y2x;

	map<gene*, gene*> s2t;
	map<gene*, gene*> t2s;

public:
	int load_genomes(const string & file1, const string & file2);
	int simulate_genomes();

	int solve();
	int solve(const string & file1, const string & file2);
	int evaluate();
	int write_mapping(const map<gene*, gene*> & m, const string & file);
};

#endif
