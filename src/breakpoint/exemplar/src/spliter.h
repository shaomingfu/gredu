#ifndef _SPLITER__H__
#define _SPLITER__H__

#include "pbase.h"
#include <vector>
#include <set>
#include <map>

using namespace std;

class spliter : public pbase
{
public:
	spliter(genome * g1, genome * g2);
	virtual ~spliter();

	vector<candidate> cps;

public:
	int split();
	int add();
	int add(const candidate & y);
	vector<candidate> compare(const candidate & x, const candidate & y) const;
	int print() const;
};

#endif
