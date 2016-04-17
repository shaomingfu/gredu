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
	: conf(_conf)
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
	}

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

int genome::write(const string & file, const string & prefix)
{
	ofstream fout(file.c_str());
	int x = 0;
	char gname[1024];
	char cname[1024];
	for(int i = 0; i < chrms.size(); i++)
	{
		if(chrms.at(i)->type == LINEAR)
		{
			gene * q = chrms.at(i)->p->b;
			while(q->b != NULL)
			{
				x++;
				sprintf(gname, "%s%d", prefix.c_str(), x);
				sprintf(cname, "chr%d", i + 1);
				fout<<gname<<" "<<q->x<<" "<<cname<<" "<<chrms[i]->type<<endl;
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
				sprintf(cname, "chr%d", i + 1);
				fout<<gname<<" "<<q->x<<" "<<cname<<" "<<chrms[i]->type<<endl;
				q = q->b;
				if(q == chrms.at(i)->p) break;
			}
		}
	}
	fout.close();
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

map<gene*, int> genome::build_gene_indices()
{
	map<gene*, int> m;
	int x = 0;
	for(int i = 0; i < chrms.size(); i++)
	{
		gene * q = chrms.at(i)->p;
		while(q != NULL)
		{
			m.insert(pair<gene*, int>(q, x));
			q = q->b;
			x++;
			if(q == chrms.at(i)->p) break;
		}
	}
	return m;
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

int genome::build_adjacencies(vector<adjacency> & s) const
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

int genome::build_adjacencies(vector<adjacency> & s, const set<gene*> & st) const
{
	s.clear();

	for(int i = 0; i < chrms.size(); i++)
	{
		chrm * ch = chrms.at(i);
		if(ch->type == LINEAR)
		{
			gene * p = ch->p;
			gene * q = p->b;
			while(q != NULL)
			{
				if(st.find(q) != st.end() || q->x == 0)
				{
					s.push_back(adjacency(PG(p, q)));
					p = q;
					q = p->b;
				}
				else
				{
					q = q->b;
				}
			}
		}
		else
		{
			assert(ch->type == CIRCULAR);
			if(ch->p == NULL) continue;

			gene * p = ch->p;
			while(true)
			{
				if(st.find(p) != st.end()) break;
				else p = p->b;
				if(p == ch->p) break;
			}

			if(st.find(p) == st.end()) continue;

			gene * a = p;
			gene * q = p->b;
			while(q != ch->p)
			{
				if(st.find(q) != st.end())
				{
					s.push_back(adjacency(PG(p, q)));
					p = q;
					q = p->b;
				}
				else
				{
					q = q->b;
				}
			}

			s.push_back(adjacency(PG(p, a)));
		}
	}

	return 0;
}


int genome::size()
{
	int s = 0;
	for(int i = 0; i < chrms.size(); i++)
	{
		s += chrms[i]->size();
	}
	return s;
}

int genome::sadist(const genome & gm, const map<gene*, gene*> & mm)
{
	map<gene*, gene*> m = mm;

	set<gene*> s1;
	set<gene*> s2;

	map<gene*, gene*>::iterator it;
	for(it = m.begin(); it != m.end(); it++)
	{
		s1.insert(it->first);
		s2.insert(it->second);
	}

	vector<adjacency> v1;
	vector<adjacency> v2;

	map<gene*, int> h;
	map<gene*, int> t;

	build_adjacencies(v1, s1);
	gm.build_adjacencies(v2, s2);
	
	for(int i = 0; i < v2.size(); i++)
	{
		adjacency a = v2[i];
		if(a.e1.b == true) h.insert(pair<gene*, int>(a.e1.g, i));
		else t.insert(pair<gene*, int>(a.e1.g, i));

		if(a.e2.b == true) h.insert(pair<gene*, int>(a.e2.g, i));
		else t.insert(pair<gene*, int>(a.e2.g, i));
	}

	int sa = 0;

	for(int i = 0; i < v1.size(); i++)
	{
		adjacency a = v1[i];
		if(a.e1.g->x != 0)
		{
			assert(m.find(a.e1.g) != m.end());
			gene * g = m[a.e1.g];
			int j = a.e1.b ? h[g] : t[g];
			if(a.weak_compare(v2[j]) >= 2) sa++;
		}
		else if(a.e2.g->x != 0)
		{
			assert(m.find(a.e2.g) != m.end());
			gene * g = m[a.e2.g];
			int j = a.e2.b ? h[g] : t[g];
			if(a.weak_compare(v2[j]) >= 2) sa++;
		}
	}

	return sa;
}

int genome::printp(const string & p) const
{
	for(int i = 0; i < chrms.size(); i++)
	{
		chrm * ch = chrms.at(i);
		if(ch->type != LINEAR) continue;
		printf("%s", p.c_str());
		gene * q = ch->p->b;
		while(q->x != 0)
		{
			if(q->b->x == 0) printf("%d", q->x);
			else printf("%d ", q->x);
			q = q->b;
		}
		printf("\n");
	}
	return 0;
}


int genome::printi(map<gene*, int> & m) const
{
	for(int i = 0; i < chrms.size(); i++)
	{
		chrm * ch = chrms.at(i);
		if(ch->type == LINEAR)
		{
			printf("[linear]      ");
			gene * q = ch->p;
			bool b = false;
			while(q != NULL)
			{
				int x = m[q];
				printf("%d(%d) ", x, q->x);
				q = q->b;
			}
		}
		else
		{
			assert(ch->type == CIRCULAR);
			if(ch->p == NULL) continue;

			printf("[circular]    ");
			gene * q = ch->p;
			bool b = false;
			while(true)
			{
				int x = m[q];
				printf("%d(%d) ", x, q->x);
				q = q->b;
			}
		}
		printf("\n");
	}
	return 0;
}


int genome::printc(int cutoff) const
{
	for(int i = 0; i < chrms.size(); i++)
	{
		chrm * ch = chrms.at(i);
		if(ch->type == LINEAR)
		{
			printf("[linear]      ");
			gene * q = ch->p;
			bool b = false;
			while(q != NULL)
			{
				if((int)fabs(q->x) <= cutoff || q->x == 0) 
				{
					printf("%d ", q->x);
					b = false;
				}
				else 
				{
					if(b == false) printf("* ");
					b = true;
				}
				q = q->b;
			}
		}
		else
		{
			assert(ch->type == CIRCULAR);
			if(ch->p == NULL) continue;

			printf("[circular]    ");
			gene * q = ch->p;
			bool b = false;
			while(true)
			{
				if((int)fabs(q->x) <= cutoff || q->x == 0) 
				{
					printf("%d ", q->x);
					b = false;
				}
				else 
				{
					if(b == false) printf("* ");
					b = true;
				}
				if(q->b == ch->p) break;
				q = q->b;
			}
		}
		printf("\n");
	}
	return 0;
}


int genome::prints(const set<gene*> & s) const
{
	for(int i = 0; i < chrms.size(); i++)
	{
		chrm * ch = chrms.at(i);
		if(ch->type == LINEAR)
		{
			printf("[linear]      ");
			gene * q = ch->p;
			while(q != NULL)
			{
				if(s.find(q) != s.end() || q->x == 0) printf("%d ", q->x);
				q = q->b;
			}
		}
		else
		{
			assert(ch->type == CIRCULAR);
			if(ch->p == NULL) continue;

			printf("[circular]    ");
			gene * q = ch->p;
			while(true)
			{
				if(s.find(q) != s.end()) printf("%d ", q->x);
				if(q->b == ch->p) break;
				q = q->b;
			}
		}
		printf("\n");
	}
	return 0;
}

bool genome::remove_tandem_linear_triple(chrm * ch)
{
	assert(ch->type == LINEAR);
	gene * x = ch->advance(ch->p, 3);
	bool result = false;
	while(x != NULL)
	{
		assert(x->a != NULL);
		assert(x->a->a != NULL);
		if(x->x + x->a->x == 0 && x->x == x->a->a->x)
		{
			operation * op = new deletion(x->a->a, x->a->a, true);
			do_deletion(op);
			delete op;
			result = true;
		}
		x = x->b;
	}
	return result;
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
		if(x->x == y->x) s++;
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
			//bool b2 = remove_tandem_linear(ch, 2);
			//bool b3 = remove_tandem_linear(ch, 3);
			//bool b4 = remove_tandem_linear(ch, 4);
			//bool b5 = remove_tandem_linear(ch, 5);
			bool b6 = remove_tandem_linear_triple(ch);
			if(!b1 && !b6) break;
		}
	}
	return 0;
}
