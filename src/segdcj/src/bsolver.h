#ifndef __BSOLVER_H__
#define __BSOLVER_H__

#include "asolver.h"
#include "shseg.h"

class bsolver : public asolver
{
public:
	bsolver(config * _conf, genome * _gm1, genome * _gm2);
	virtual ~bsolver();

public:
	vector<shseg> vss;	// shared segments
	set<int> sdup;		// list of possibly duplicated genes
	ugraph mr;			// graph for checking alternating path
	vector<int> mvim;	// the index map of the matching graph
	vector<int> mate;	// initial matching

public:
	// build possibly duplicated genes
	int build_possibly_duplicated_genes();

	// matching graph and checking alternating path
	int build_matching_graph();
	int build_initial_matching();
	bool check_alternating_path(int s, int t);

	// build shared segments
	int build_shared_segments();
	int add_shared_segment(int x, int y);
	int extend_forward(int x, int y, vector<int> & xv, vector<int> & yv, set<int> & sx, set<int> & sy);
	int extend_backward(int x, int y, vector<int> & xv, vector<int> & yv, set<int> & sx, set<int> & sy);

	// make segments
	int make_shared_segments();
	int make_shared_segment(shseg & ss);

	// prepare segments
	int prepare_shared_segments();
	int prepare_shared_segment(shseg & ss);

	// verity segments
	int verify_shared_segments();
	bool verify_shared_segment(shseg & ss);

	// fix segment
	int fix_shared_segment(const shseg & ss);

	// print
	int print_shared_segments();
	int assert_mate();
	int statistic();

	// draw tex
	int draw_adjacency_graph(const string & file, shseg & ss);
};

#endif
