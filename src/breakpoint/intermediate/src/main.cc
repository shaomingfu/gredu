#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cstdio>

#include "trier.h"
#include "mygraph.h"

using namespace std;

int main(int argc, const char ** argv)
{
	if(argc != 4)
	{
		cout<<"usage: "<<argv[0]<<" <genome1> <genome2> [ilp-timelimit]"<<endl;
		return 0;
	}

	srand(time(0));

	trier tr(atof(argv[3]));
	tr.solve(argv[1], argv[2]);

	return 0;
}
