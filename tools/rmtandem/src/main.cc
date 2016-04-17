#include <iostream>
#include <cstdlib>
#include <cstdio>

#include "genome.h"

using namespace std;

int main(int argc, const char ** argv)
{
	if(argc != 3)
	{
		cout<<"usage: "<<argv[0]<<" <input-genome> <output-genome>"<<endl;
		return 0;
	}

	genome g;
	g.load(argv[1]);
	g.remove_tandem();
	g.write(argv[2], "");

	return 0;
}
