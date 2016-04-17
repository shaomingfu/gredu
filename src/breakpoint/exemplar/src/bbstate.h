#ifndef __BBSTATE_H__
#define __BBSTATE_H__

#include "gene.h"
#include <set>
#include <vector>
#include <map>

using namespace std;

class bbstate
{
public:
	bbstate(const vector<gene*> & _ig);

public:
	const vector<gene*> & ig;

	set<int> sx;					// current genes (id) in genome1
	set<int> sy;					// current genes (id) in genome2
	vector<int> fx;					// -1 if this family is not in; otherwise it is the gene (id)
	vector<int> fy;					// -1 if this family is not in; otherwise it is the gene (id)
	int distance;					// number of breakpoint in this current partial solution

	vector<bool> sab;				// false: the shared adjacency is destroyed; true: still usable
	vector<bool> cdb;				// false: the candidate is fixed; true: the candidate is still usable
	map< int, set<int> > cdm;		// key: shadj index, value: affected candidates

public:
	int affect_candidate(int sa, int cd);
	int add_pair(const PI & p);
	int add_pair(const PI & p, int inc);
	int remove_pair(const PI & p);
	int remove_pair(const PI & p, int dec);
	int inc_distance(const PI & p) const;
	bool adjacent(int xl, int xu, int yl, int yu) const;
	int print() const;
};

int find_less(const set<int> & s, int x);
int find_less_equal(const set<int> & s, int x);
int find_greater(const set<int> & s, int x);
int find_greater_equal(const set<int> & s, int x);
bool neighbor(const set<int> & s, int x, int y);

#endif
