#include "config.h"
#include <cstdlib>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cstring>

using namespace std;

config::config()
{
	ilp_focus = 1;
	concurrent = 1;
	algo = 1;
	num_exemplars = 1;
	max_combination = 100;
	heuristic = false;
}

config::config(const char * conf_file)
{
	load_config(conf_file);
}

int config::load_config(const char * conf_file)
{
	ifstream fin(conf_file);
	if(fin.fail())
	{
		cout<<"open config file error "<<conf_file<<endl;
		return -1;
	}

	char buf[1024];
	char key[1024];
	char value[1024];
	
	while(fin.getline(buf, 1024, '\n'))
	{
		stringstream sstr(buf);
		sstr>>key>>value;

		if(strcmp(key, "alphabet_size")==0)
		{
			alphabet_size = atoi(value);
		}
		else if(strcmp(key, "min_duplication_length")==0)
		{
			min_duplication_length = atoi(value);
		}
		else if(strcmp(key, "max_duplication_length")==0)
		{
			max_duplication_length = atoi(value);
		}
		else if(strcmp(key, "simulation_length")==0)
		{
			simulation_length = atoi(value);
		}
		else if(strcmp(key, "random_duplication_num")==0)
		{
			random_duplication_num = atoi(value);
		}
		else if(strcmp(key, "prob_inversion")==0)
		{
			prob_inversion = atof(value);
		}
		else if(strcmp(key, "prob_duplication")==0)
		{
			prob_duplication = atof(value);
		}
		else if(strcmp(key, "ilp_timelimit")==0)
		{
			ilp_timelimit = atof(value);
		}
		else if(strcmp(key, "concurrent")==0)
		{
			concurrent = atoi(value);
		}
		else if(strcmp(key, "num_exemplars")==0)
		{
			num_exemplars = atoi(value);
		}
		else if(strcmp(key, "algo")==0)
		{
			algo = atoi(value);
		}
		else if(strcmp(key, "ilp_focus")==0)
		{
			ilp_focus = atoi(value);
		}
		else if(strcmp(key, "max_combination")==0)
		{
			max_combination = atoi(value);
		}
		else if(strcmp(key, "heuristic")==0)
		{
			if(strcmp(value, "true") == 0) heuristic = true;
			else heuristic = false;
		}
	}
	return 0;
}

