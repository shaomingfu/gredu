#ifndef __CONFIG_H__
#define __CONFIG_H__

#include <string>
using namespace std;

class config
{
public:
	config();
	config(const char * conf_file);
	int load_config(const char * conf_file);

public:
	int alphabet_size;

	int max_duplication_length;
	int min_duplication_length;
	int random_duplication_num;
	int simulation_length;

	double prob_inversion;
	double prob_duplication;

	double ilp_timelimit;
	int ilp_focus;
	int concurrent;

	bool remove_redundant;
	bool use_psolver;
	int max_combination;
	bool heuristic;
	int algo;
};

#endif
