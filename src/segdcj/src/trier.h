#ifndef __TRIER_H__
#define __TRIER_H__

#include "genome.h"
#include "config.h"

typedef map<gene*, gene*> MPG;

class trier
{
public:
	trier(const char * file);
	trier(double limit);
	virtual ~trier();
	
public:
	config * conf;
	genome * gm1;
	genome * gm2;

public:
	int load_genomes(const string & file1, const string & file2);
	int simulate_genomes();

	int solve();
	int solve(const string & file1, const string & file2);
	int write_mapping(const MPG & x2y, const string & file);
};

#endif
