#ifndef __CONFIG_H__
#define __CONFIG_H__

#include <string>
using namespace std;

class config
{
public:
	config(double time);

public:
	int alphabet_size;
	double ilp_timelimit;
};

#endif
