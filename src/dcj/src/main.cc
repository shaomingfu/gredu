#include <iostream>
#include <cstdlib>
#include <cstdio>

#include "gredo.h"

using namespace std;

int main(int argc, const char ** argv)
{
	if(argc != 4 && argc != 3)
	{
		cout<<"usage: "<<argv[0]<<" <genome1> <genome2> [ILP timelimit (default=7200 seconds)]"<<endl;
		return 0;
	}

	double time = 7200;
	if(argc == 4) time = atoi(argv[3]);

	gredo g(time);
	g.solve(argv[1], argv[2]);

	return 0;
}
