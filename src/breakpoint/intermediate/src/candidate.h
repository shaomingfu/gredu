#ifndef __CANDIDATE_H__
#define __CANDIDATE_H__

#include "gene.h"
#include <vector>
#include <set>

using namespace std;

class candidate
{
public:
	candidate();
	candidate(bool b, const vector<gene*> & xv, const vector<gene*> & yv);

public:
	bool d;
	vector<gene*> xv;
	vector<gene*> yv;
	//set<int> fs;

public:
	int print() const;
	int size() const;
	bool operator< (const candidate & c) const;
	
	//bool conflict(const candidate & c) const;
	//int split(int p, candidate & x, candidate & y) const;
};

#endif
