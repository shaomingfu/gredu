#ifndef __DECOMPOSER_H__
#define __DECOMPOSER_H__

#include "pbase.h"
#include "mygraph.h"

using namespace std;

typedef pair<int, int> PI;
typedef map<int, int> MPI;

class decomposer : public pbase
{
public:
	decomposer(genome * g1, genome * g2);
	virtual ~decomposer();

public:
	ugraph gr;
	MPI f2g;
	MPI g2f;

public:
	int solve();
	int build();
	int print_neato(const string & file);
	vector<int> collect(gene * x, gene * y);
};

#endif
