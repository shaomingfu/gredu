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

bool genome::remove_tandem_no_sign(chrm * ch, int p)
{
	assert(ch->type == LINEAR);
	gene * x = ch->p->b;
	gene * y = ch->advance(x, p);
	gene * z = NULL;
	bool result = false;
	int s = 0;
	while(x != NULL && y != NULL)
	{
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
			bool b0 = remove_tandem_no_sign(ch, 1);
			//bool b1 = remove_tandem_linear(ch, 1);
			//bool b2 = remove_tandem_linear(ch, 2);
			//bool b3 = remove_tandem_linear(ch, 3);
			//bool b4 = remove_tandem_linear(ch, 4);
			//bool b5 = remove_tandem_linear(ch, 5);
			//bool b6 = remove_tandem_linear_triple(ch);
			if(!b0) break;
		}
	}
	return 0;
}
