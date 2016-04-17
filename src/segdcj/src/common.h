#ifndef __COMMON_H__
#define __COMMON_H__

#include <vector>
#include <cassert>

using namespace std;

const static int OTHER = -1;
const static int EQUAL = 1;
const static int SUBSET = 2;
const static int REVERSE_EQUAL = 3;
const static int REVERSE_SUBSET = 4;

const static int EMPTY = -1;

const static int GENE = 0;
const static int ADJACENCY = 1;

const static int SUPER = -3;

const static int START = -1;
const static int END = -2;

const static int ACTIVE = 0;
const static int NONACTIVE = 1;
const static int NONACTIVE_SOURCE = 2;
const static int NONACTIVE_TARGET = 3;

const static int LINEAR = 1;
const static int CIRCULAR = 2;

const static int INVERSION = 1;
const static int DUPLICATION = 2;
const static int LOSS = 3;
const static int INSERTION = 4;
const static int DELETION = 5;
const static int IDUPLICATION = 6;
const static int CIRCULARIZATION = 7;
const static int EXCHANGE = 8;
const static int DCJ = 9;
const static int REVERSION = 10;
const static int LINEARIZATION = 11;
const static int CDUPLICATION = 12;

const static double ZERO = 1e-10;
const static double MAX_NUM = 1e30;
const static double MIN_NUM = -1e30;

double log_add(double x, double y);

int random(const vector<double> & x);

int invert(vector<int> & v);

template<typename T>
int combine(vector<T> & v1, const vector<T> & v2)
{
	for(int i = 0; i < v2.size(); i++)
	{
		v1.push_back(v2.at(i));
	}
	return 0;
}

template<typename T>
int shortcut(vector<T> & v, int x, int y)
{
	assert(x >= 0 && x < v.size());
	assert(y >= 0 && y < v.size());
	assert(x <= y);

	for(int i = y + 1; i < v.size(); i++)
	{
		int k = i + x - y - 1;
		v.at(k) = v.at(i);
	}

	v.resize(v.size() + x - y - 1);
	return 0;
}

#endif
