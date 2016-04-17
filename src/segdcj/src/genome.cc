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
	: conf(_conf), genome_base()
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


int genome::count_inversion()
{
	int sum = 0;
	for(int i = 0; i < chrms.size(); i++)
	{
		sum += chrms.at(i)->count_inversion();
	}
	return sum;
}

int genome::count_duplication()
{
	int sum = 0;
	for(int i = 0; i < chrms.size(); i++)
	{
		sum += chrms.at(i)->count_duplication(conf->max_duplication_length);
	}
	return sum;
}


int genome::count_loss()
{
	int sum = 0;
	for(int i = 0; i < chrms.size(); i++)
	{
		sum += chrms.at(i)->count_loss();
	}
	return sum;
}

int genome::enumerate_operations(vector<operation*> & operations) const
{
	assert(operations.size()==0);

	for(int i = 0; i < chrms.size(); i++)
	{
		chrms.at(i)->enumerate_reverse_duplication(operations, conf->max_duplication_length);
		chrms.at(i)->enumerate_reverse_iduplication(operations, conf->max_duplication_length);
		for(int j = i + 1; j < chrms.size(); j++)
		{
			chrms.at(i)->enumerate_reverse_duplication(chrms[j], operations, conf->max_duplication_length);
			chrms.at(i)->enumerate_reverse_iduplication(chrms[j], operations, conf->max_duplication_length);
		}
	}

	// TODO, ignore templates of the duplications
	set<PG> s;
	vector<operation*> ops;
	for(int i = 0; i < operations.size(); i++)
	{
		operation * op = operations[i];
		gene * head = NULL;
		gene * tail = NULL;
		if(op->type == DUPLICATION)
		{
			duplication * opx = dynamic_cast<duplication*>(op);
			head = opx->head;
			tail = opx->tail;
		}
		else if(op->type == IDUPLICATION)
		{
			iduplication * opx = dynamic_cast<iduplication*>(op);
			head = opx->head;
			tail = opx->tail;
		}

		PG p(head, tail);
		if(s.find(p) == s.end())
		{
			ops.push_back(op);
			s.insert(p);
		}
		else
		{
			delete op;
		}
	}

	operations.clear();
	operations = ops;

	return 0;
}

int genome::random_duplicate()
{
	int c = rand() % chrms.size();

	int max_length = conf->max_duplication_length;
	int min_length = conf->min_duplication_length;

	if(max_length > chrms.at(c)->size())
	{
		max_length = chrms.at(c)->size();
	}

	int len = min_length + (rand() % (max_length - min_length + 1));
	//int len = rand() % conf->max_duplication_length + 1;

	//printf("duplication length = %4d\n", len);

	if(chrms.at(c)->type==LINEAR)
	{
		int start = rand() % ( chrms.at(c)->size() - len + 1 );
		int dest = rand() % ( chrms.at(c)->size() + 1 );

		gene * gstart = chrm::advance(chrms.at(c)->p->b, start);
		gene * gend = chrm::advance(gstart, len - 1);
		gene * gdest = chrm::advance(chrms.at(c)->p->b, dest);

		operation * op = new duplication(gdest, gstart, gend);

		operate(op);
		delete op;
	}
	else
	{
		int start = rand() % (chrms.at(c)->size());
		int dest = rand() % (chrms.at(c)->size());

		gene * gstart = chrm::advance(chrms.at(c)->p, start);
		gene * gend = chrm::advance(gstart, len - 1);
		gene * gdest = chrm::advance(chrms.at(c)->p, dest);

		operation * op = new duplication(gdest, gstart, gend);

		operate(op);
		delete op;
	}

	return 0;
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

map<gene*, int> genome::build_gene_indices()
{
	map<gene*, int> gi;

	int index = 0;
	for(int i = 0; i < chrms.size(); i++)
	{
		chrm * ch = chrms.at(i);
		if(ch->type == LINEAR)
		{
			gene * q = ch->p;
			while(true)
			{
				gi.insert(pair<gene*, int>(q, index));
				index++;
				if(q->b == NULL) break;
				q = q->b;
			}
		}
		else
		{
			gene * q = ch->p;
			while(true)
			{
				assert(q->x != 0);
				gi.insert(pair<gene*, int>(q, index));
				q = q->b;
				if(q == ch->p) break;
			}
		}
	}
	return gi;
}

int genome::build_gene_indices(map<gene*, int> & gi, map<int, gene*> & ig)
{
	gi.clear();
	ig.clear();

	int index = 0;
	for(int i = 0; i < chrms.size(); i++)
	{
		chrm * ch = chrms.at(i);
		if(ch->type == LINEAR)
		{
			gene * q = ch->p;
			while(true)
			{
				gi.insert(pair<gene*, int>(q, index));
				ig.insert(pair<int, gene*>(index, q));
				index++;
				if(q->b == NULL) break;
				q = q->b;
			}
		}
		else
		{
			gene * q = ch->p;
			while(true)
			{
				assert(q->x != 0);
				gi.insert(pair<gene*, int>(q, index));
				ig.insert(pair<int, gene*>(index, q));
				index++;
				q = q->b;
				if(q == ch->p) break;
			}
		}
	}
	return 0;
}

vector<int> genome::build_gene_copy()
{
	vector<int> list;
	list.assign(conf->alphabet_size + 1, 0);

	for(int i = 0; i < chrms.size(); i++)
	{
		chrm * ch = chrms.at(i);
		if(ch->type == LINEAR)
		{
			gene * q = ch->p;
			while(q != NULL)
			{
				int f = (int)fabs(q->x);
				if(f <= conf->alphabet_size) list.at(f)++;
				q = q->b;
			}
		}
		else
		{
			gene * q = ch->p;
			while(true)
			{
				assert(q->x != 0);
				int f = (int)fabs(q->x);
				if(f <= conf->alphabet_size) list.at(f)++;

				list.at( abs(q->x) )++;
				q = q->b;
				if(q == ch->p) break;
			}
		}
	}
	return list;
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
			int f = (int)fabs(q->x);
			if(f <= conf->alphabet_size) list.at(f).push_back(q);
			q = q->b;
			if(q == chrms.at(i)->p) break;
		}
	}
}

int genome::build_gene_map(vector< vector<gene*> > & list, gene * head, gene * tail)
{
	assert(head->ch == tail->ch);

	list.clear();

	list.resize(conf->alphabet_size + 1);

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

int genome::write(const string & file)
{
	ofstream fout(file.c_str());
	char cname[1024];
	for(int i = 0; i < chrms.size(); i++)
	{
		sprintf(cname, "chr%d", i + 1);
		if(chrms.at(i)->type == LINEAR)
		{
			gene * q = chrms.at(i)->p->b;
			while(q->b != NULL)
			{
				fout<<q->s.c_str()<<" "<<q->x<<" "<<cname<<" "<<chrms[i]->type<<endl;
				q = q->b;
			}
		}
		else if(chrms.at(i)->type == CIRCULAR)
		{
			gene * q = chrms.at(i)->p;
			while(true)
			{
				fout<<q->s.c_str()<<" "<<q->x<<" "<<cname<<" "<<chrms[i]->type<<endl;
				q = q->b;
				if(q == chrms.at(i)->p) break;
			}
		}
	}
	fout.close();
	return 0;
}

int genome::set_names(const string & prefix)
{
	int x = 0;
	char gname[1024];
	for(int i = 0; i < chrms.size(); i++)
	{
		if(chrms.at(i)->type == LINEAR)
		{
			gene * q = chrms.at(i)->p->b;
			while(q->b != NULL)
			{
				x++;
				sprintf(gname, "%s%d", prefix.c_str(), x);
				q->s = gname;
				q = q->b;
			}
		}
		else if(chrms.at(i)->type == CIRCULAR)
		{
			gene * q = chrms.at(i)->p;
			while(true)
			{
				x++;
				sprintf(gname, "%s%d", prefix.c_str(), x);
				q->s = gname;
				q = q->b;
				if(q == chrms.at(i)->p) break;
			}
		}
	}
	return 0;
}

set<gene*> genome::merge_tandem()
{
	set<gene*> s;
	for(int i = 0; i < chrms.size(); i++)
	{
		chrm * ch = chrms.at(i);
		merge_tandem(ch, s);
	}
	return s;
}

int genome::merge_tandem(chrm * ch, set<gene*> & s)
{
	gene * x = ch->p;
	gene * z = NULL;
	while(x != NULL && x->b != NULL)
	{
		if((x->x == x->b->x || x->x + x->b->x == 0) && x->x != 0)
		{
			if(z == NULL) z = x;
			x = x->b;
		}
		else if(z != NULL)
		{
			x = x->b;
			s.insert(z);
			operation * op = new deletion(z->b, x->a, true);
			do_deletion(op);
			delete op;
			z = NULL;
		}
		if(x == ch->p) break;	
	}
	return 0;
}

bool genome::remove_tandem_linear(chrm * ch, int p)
{
	assert(ch->type == LINEAR);
	gene * x = ch->p->b;
	gene * y = ch->advance(x, p);
	gene * z = NULL;
	bool result = false;
	int s = 0;
	while(x != NULL && y != NULL)
	{
		// TODO
		//if(x->x == y->x) s++;
		if(x->x == y->x || x->x + y->x == 0) s++;
		else s = 0;
		if(s == 1) z = x;
		x = x->b;
		y = y->b;
		if(s < p) continue;
		operation * op = new deletion(z, x->a, true);
		do_deletion(op);
		delete op;
		result = true;
		s = 0;
	}
	return result;
}

int genome::remove_tandem()
{
	for(int i = 0; i < chrms.size(); i++)
	{
		chrm * ch = chrms.at(i);
		if(ch->type == CIRCULAR) continue;
		while(true)
		{
			bool b1 = remove_tandem_linear(ch, 1);
			bool b2 = remove_tandem_linear(ch, 2);
			bool b3 = remove_tandem_linear(ch, 3);
			bool b4 = remove_tandem_linear(ch, 4);
			bool b5 = remove_tandem_linear(ch, 5);
			if(!b1 && !b2 && !b3 && !b4 && !b5) break;
		}
	}
	return 0;
}

int genome::statistic()
{
	vector<int> v = build_gene_copy();
	int empty = 0;
	int single = 0;
	int multi = 0;
	int total = 0;

	for(int i = 0; i < v.size(); i++)
	{
		total += v[i];
		if(v[i] == 0) empty++;
		else if(v[i] == 1) single++;
		else multi++;
	}
	printf("total (%6d/%5d) families/genes, %6d empty ones, %6d singletons, (%6d,%6d) multiple-gene families/genes\n",
			(int)v.size(), total, empty, single, multi, total - single);
	return 0;
}
