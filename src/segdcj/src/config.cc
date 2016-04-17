#include "config.h"
#include <cstdlib>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cstring>

using namespace std;

config::config(double limit)
{
	ilp_timelimit = limit;
	ilp_focus = 1;
	concurrent = 1;
	sd_weight = 0.75;
	heuristics = 0.5;
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
		else if(strcmp(key, "random_duplication_num")==0)
		{
			random_duplication_num = atoi(value);
		}
		else if(strcmp(key, "inversion1")==0)
		{
			inversion1 = atoi(value);
		}
		else if(strcmp(key, "duplication1")==0)
		{
			duplication1 = atoi(value);
		}
		else if(strcmp(key, "inversion2")==0)
		{
			inversion2 = atoi(value);
		}
		else if(strcmp(key, "duplication2")==0)
		{
			duplication2 = atoi(value);
		}
		else if(strcmp(key, "ilp_timelimit")==0)
		{
			ilp_timelimit = atof(value);
		}
		else if(strcmp(key, "heuristics")==0)
		{
			heuristics = atof(value);
		}
		else if(strcmp(key, "concurrent")==0)
		{
			concurrent = atoi(value);
		}
		else if(strcmp(key, "ilp_focus")==0)
		{
			ilp_focus = atoi(value);
		}
		else if(strcmp(key, "sd_weight")==0)
		{
			sd_weight = atof(value);
		}
		else if(strcmp(key, "algo")==0)
		{
			algo = atoi(value);
		}
	}
	return 0;
}

