#ifndef __SHSEG_H__
#define __SHSEG_H__

#include "mygraph.h"

#include <boost/graph/connected_components.hpp>
#include <vector>
#include <set>

using namespace std;

typedef pair<int, int> PI;
typedef vector< vector<int> > VVI;

class shseg
{
public:
	shseg();
	shseg(const vector<int> & x, const vector<int> & y);
	shseg(const vector<PI> & v);

public:
	vector<int> xv;		// x part of the segment
	vector<int> yv;		// y part of the segment

	set<int> sa;		// all extremities connecting to the segment
	set<int> sb;		// adjacent extremities to sa
	set<PI> ea;			// edges between the extremities in sa
	set<PI> eb;			// possible edges between extremities in sb
	PI xb;				// boundary of xv
	PI yb;				// boundary of yv

	map<int, int> ei;	// indices for extremities (in sa and sb) 
	map<int, int> ie;	// indices for extremities (in sa and sb) 
	ugraph ar;			// adjacency graph
	vector<PI> ge;		// list of genes (head and tail)
	vector<int> gf;		// gene families
	vector<int> adj;	// adjacent extremties
	VVI vdup;			// duplicons

public:
	bool operator< (const shseg & ss) const;
	int size() const;
	int calc_lower_bound();
	int statistic() const;
	int print() const;
	bool contain(int x, int y) const;

	map<int, int> sort_vertices(int xsize, bool b);
	int draw(const string & file, int xsize);
};

#endif
