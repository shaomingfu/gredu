#include "config.h"
#include <cstdlib>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cstring>

using namespace std;

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
		else if(strcmp(key, "max_duplication_length")==0)
		{
			max_duplication_length = atoi(value);
		}
		else if(strcmp(key, "max_sampling_length")==0)
		{
			max_sampling_length = atoi(value);
		}
		else if(strcmp(key, "simulation_length")==0)
		{
			simulation_length = atoi(value);
		}
		else if(strcmp(key, "random_duplication_num")==0)
		{
			random_duplication_num = atoi(value);
		}
		else if(strcmp(key, "sample_path_num")==0)
		{
			sample_path_num = atoi(value);
		}


	}
	return 0;
}

