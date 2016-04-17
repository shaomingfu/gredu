#ifndef __GENE_H__
#define __GENE_H__

#include <string>
#include <utility>
#include <cstddef>
#include <vector>

using namespace std;

class chrm;

class gene
{
public:
	gene(int _x);
	gene(int _x, chrm * _ch);
	gene(int _x, chrm * _ch, const string & _s);
	//gene(const gene & g);

public:
	int x;
	string s;
	gene * a;
	gene * b;
	chrm * ch;

public:
	int reverse();
	static int get_gene_copy_index(gene * g);
};

typedef pair<gene*, gene*> PG;
typedef pair<int, int> PI;

vector<gene*> build_gene_list(const PG & pg);
int distance(const PG & pg);

#endif
