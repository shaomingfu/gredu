#ifndef __SIMULATOR_H__
#define __SIMULATOR_H__

#include "genome.h"
#include "config.h"

#include <map>

using namespace std;

class simulator
{
public:
	simulator(config * _conf, genome * _target1, genome * _target2);
	~simulator();

public:
	genome source;
	genome * target1;
	genome * target2;

	config * conf;

	map<gene*, gene*> s2t;
	map<gene*, gene*> t2s;

public:
	int init_source();

	int simulate();
	int simulate1();
	int simulate2();

	operation * simulate_single(genome * t);
	operation * simulate_dcj(genome * t);
	operation * simulate_dup(genome * t);
};

#endif
