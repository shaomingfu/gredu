#ifndef __CONFIG_H__
#define __CONFIG_H__

#include <string>
using namespace std;

class config
{
public:
	config(double limit);
	config(const char * conf_file);

	int load_config(const char * conf_file);

public:
	int alphabet_size;
	int max_duplication_length;
	int min_duplication_length;
	int random_duplication_num;

	int inversion1;
	int inversion2;
	int duplication1;
	int duplication2;

	double ilp_timelimit;
	int ilp_focus;
	int concurrent;
	double heuristics;

	double sd_weight;
	int algo;
};

#endif
