#ifndef __ADJACENCY__
#define __ADJACENCY__

#include <vector>
#include <string>
#include "extremity.h"

using namespace std;

class adjacency
{
public:
	adjacency(const adjacency & a);
	adjacency();
	adjacency(const PG & pg);
	adjacency(const extremity & _e1, const extremity & _e2);

	extremity e1;
	extremity e2;

public:
	adjacency & operator=(const adjacency & a);
	bool operator==(const adjacency & a) const;
	bool operator!=(const adjacency & a) const;
	bool available() const;

	string label_str() const;
	string label_tex() const;

	int print() const;
	bool direction() const;
	int exchange();
	int weak_compare(const adjacency & a) const;
	int strong_compare(const adjacency & a) const;

	adjacency diff(const adjacency & a, const adjacency & o) const;
	extremity diff(const adjacency & o) const;
	extremity same(const adjacency & o) const;

	static int sort_adjacencies(vector<adjacency> & s);
};

vector<adjacency> build_adjacency_list(const PG & pg);

#endif
