#include "genome.h"
#include "chrm.h"
#include "linear_chrm.h"
#include "circular_chrm.h"

#include <algorithm>
#include <sstream>
#include <fstream>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <set>
#include <cstring>

using namespace std;

genome::genome(config * _conf)
	: conf(_conf), genome_base(_conf->alphabet_size)
{
}

genome& genome::operator=(const genome & g)
{
	assert(g.conf == this->conf);
	genome_base::operator=(g);
	return (*this);
}

/*
genome::genome(const genome_base & g, config * _conf)
	: genome_base(g)
{
	conf = _conf;
}
*/

genome::genome(const genome & g)
	: genome_base(g)
{
	conf = g.conf;
}

genome::~genome()
{
	// TODO
	//genome_base::~genome_base();
}

int genome::load(const string & file)
{
	int num = 0;

	ifstream fin(file.c_str());

	if(fin.fail())
	{
		printf("open file error: %s\n", file.c_str());
		return -1;
	}

	char buf[10240];
	vector<int> v;
	vector<string> s;
	string pname = "ASDFGHJKL";
	int ptype = -1;

	char gname[10240];
	int family = -1;
	char cname[10240];
	int type = -1;

	while(fin.getline(buf, 10240, '\n'))
	{
		stringstream sstr(buf);
		sstr>>gname>>family>>cname>>type;
		if(string(cname) != pname && v.size() > 0)
		{
			if(ptype == LINEAR) add_linear_chrm(v, s);
			else if(ptype == CIRCULAR) add_circular_chrm(v, s);
			else assert(false);
			v.clear();
			s.clear();
		}

		ptype = type;
		pname = cname;
		v.push_back(family);
		s.push_back(gname);
	}
	if(v.size() > 0)
	{
		if(ptype == LINEAR) add_linear_chrm(v, s);
		else if(ptype == CIRCULAR) add_circular_chrm(v, s);
		else assert(false);
	}

	fin.close();
	return 0;
}

int genome::build_gene_map(vector< vector<gene*> > & list)
{
	list.clear();

	list.resize(conf->alphabet_size + 1);

	for(int i = 0; i < chrms.size(); i++)
	{
		gene * q = chrms.at(i)->p;
		while(q != NULL)
		{
			list.at( abs(q->x) ).push_back(q);
			q = q->b;
			if(q == chrms.at(i)->p) break;
		}
	}
}

int genome::find_max_gene()
{
	int max = 0;
	for(int i = 0; i < chrms.size(); i++)
	{
		chrm * ch = chrms.at(i);
		if(ch->type == LINEAR)
		{
			gene * q = ch->p->b;
			while(true)
			{
				assert(q->x != 0);
				max = max > abs(q->x) ? max : abs(q->x);
				q = q->b;
				if(q->b == NULL) break;
			}
		}
		else
		{
			gene * q = ch->p;
			while(true)
			{
				assert(q->x != 0);
				max = max > abs(q->x) ? max : abs(q->x);
				q = q->b;
				if(q == ch->p) break;
			}
		}
	}
	return max;
}

vector<int> genome::build_gene_copy()
{
	vector<int> list;
	list.assign(conf->alphabet_size, 0);

	for(int i = 0; i < chrms.size(); i++)
	{
		chrm * ch = chrms.at(i);
		if(ch->type == LINEAR)
		{
			gene * q = ch->p->b;
			while(true)
			{
				assert(q->x != 0);
				list.at( abs(q->x) - 1 )++;
				q = q->b;
				if(q->b == NULL) break;
			}
		}
		else
		{
			gene * q = ch->p;
			while(true)
			{
				assert(q->x != 0);
				list.at( abs(q->x) - 1 )++;
				q = q->b;
				if(q == ch->p) break;
			}
		}
	}
	return list;
}

int genome::build_gene_map(vector< vector<gene*> > & list, gene * head, gene * tail)
{
	assert(head->ch == tail->ch);

	list.clear();

	list.resize(alphabet_size + 1);

	for(int i = 0; i < chrms.size(); i++)
	{
		if(chrms.at(i) == head->ch) continue;

		gene * q = chrms.at(i)->p;
		while(q != NULL)
		{
			list.at( abs(q->x) ).push_back(q);
			q = q->b;
			if(q == chrms.at(i)->p) break;
		}
	}

	if(head->ch->type == LINEAR)
	{
		gene * q = tail->b;
		while(q != NULL)
		{
			list.at( abs(q->x) ).push_back(q);
			q = q->b;
		}

		q = head->a;
		while(q != NULL)
		{
			list.at( abs(q->x) ).push_back(q);
			q = q->a;
		}
	}
	else
	{
		gene * q = tail->b;
		while(true)
		{
			if(q == head) break;
			list.at( abs(q->x) ).push_back(q);
			q = q->b;
		}
	}

	return 0;
}

int genome::build_adjacencies(vector<adjacency> & s)
{
	s.clear();

	for(int i = 0; i < chrms.size(); i++)
	{
		chrm * ch = chrms.at(i);
		if(ch->type == LINEAR)
		{
			gene * q = ch->p;
			while(true)
			{
				if(q->b == NULL)
				{
					assert(q->x == 0);
					break;
				}

				s.push_back(adjacency(PG(q, q->b)));
				q = q->b;
			}
		}
		else
		{
			assert(ch->type == CIRCULAR);
			if(ch->p == NULL) continue;

			gene * q = ch->p;
			while(true)
			{
				s.push_back(adjacency(PG(q, q->b)));
				if(q->b == ch->p) break;
				q = q->b;
			}
		}
	}

	return 0;
}

