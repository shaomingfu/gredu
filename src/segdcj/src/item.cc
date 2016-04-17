#include "item.h"
#include "common.h"

#include <cmath>
#include <cassert>
#include <cstdio>

item * copy_item(item * x)
{
	if(x->type == GENE)
	{
		gene_item * xx = dynamic_cast<gene_item*>(x);
		return new gene_item(xx->g);
	}
	else if(x->type == ADJACENCY)
	{
		adj_item * xx = dynamic_cast<adj_item*>(x);
		return new adj_item(xx->a);
	}
	else if(x->type == CDUPLICATION)
	{
		cdup_item * xx = dynamic_cast<cdup_item*>(x);
		return new cdup_item(*xx);
	}
	else if(x->type == DUPLICATION)
	{
		dup_item * xx = dynamic_cast<dup_item*>(x);
		return new dup_item(*xx);
	}
	else if(x->type == DCJ)
	{
		dcj_item * xx = dynamic_cast<dcj_item*>(x);
		return new dcj_item(*xx);
	}
	else assert(1 == 0);
}

gene_item::gene_item(gene * _g)
{
	type = GENE;
	g = _g;
	valid = true;
	make_name();
}

int gene_item::make_name()
{
	//name = "G";
	char buf[1024];
	sprintf(buf, "%d.%d", (int)(fabs(g->x)), gene::get_gene_copy_index(g));
	name.append("(");
	name.append(buf);
	name.append(")");
	return 0;
}

int gene_item::find_adjacency(const adjacency & x)
{
	return -1;
}

int gene_item::find_gene(gene * _g)
{
	if(g == _g) return 0;
	else return -1;
}

int gene_item::operate(genome_base & gm) const
{
	return 0;
}

int gene_item::reverse_operate(genome_base & gm) const
{
	return 0;
}

int gene_item::print() const
{
	printf("gene item: %d\n", g->x);
	return 0;
}

int gene_item::check_in_active(gene * g)
{
	return NONACTIVE;
}

int gene_item::check_out_active(gene * g)
{
	return NONACTIVE;
}

int gene_item::check_in_active(const adjacency & a)
{
	return NONACTIVE;
}

int gene_item::check_out_active(const adjacency & a)
{
	return NONACTIVE;
}

adj_item::adj_item(const adjacency & _a)
{
	type = ADJACENCY;
	a = _a;
	valid = true;
	make_name();
}

int adj_item::make_name()
{
	//name = "A";
	name = "";
	name.append(a.label_str());
	return 0;
}

int adj_item::find_adjacency(const adjacency & x)
{
	if(x == a) return 0;
	return -1;
}

int adj_item::find_gene(gene * _g)
{
	return -1;
}

int adj_item::operate(genome_base & gm) const
{
	return 0;
}

int adj_item::reverse_operate(genome_base & gm) const
{
	return 0;
}

int adj_item::print() const
{
	printf("adj item: ");
	a.print();
	return 0;
}

int adj_item::check_in_active(gene * g)
{
	return NONACTIVE;
}

int adj_item::check_out_active(gene * g)
{
	return NONACTIVE;
}

int adj_item::check_in_active(const adjacency & a)
{
	return NONACTIVE;
}

int adj_item::check_out_active(const adjacency & a)
{
	return NONACTIVE;
}

extremity cdup_item::get_target_boundary_extremity(const extremity & e)
{
	// TODO, assume that there is one-to-one correspondence between sg and tg
	for(int i = 0; i < sg.size(); i++)
	{
		assert((int)fabs(sg.at(i)->x) == (int)fabs(tg.at(i)->x));
		if(sg.at(i) == e.g) return extremity(tg.at(i), e.b);
	}
	assert(1 == 0);
}

cdup_item::cdup_item(const vector<adjacency> & _s, const vector<adjacency> & _t, const adjacency & _c, const vector<gene*>& _sg, const vector<gene*> & _tg)
{
	s = _s;
	t = _t;
	c = _c;
	sg = _sg;
	tg = _tg;
	type = CDUPLICATION;
	valid = true;
	make_name();
}
cdup_item::cdup_item(const cdup_item & ci)
{
	s = ci.s;
	t = ci.t;
	c = ci.c;
	sg = ci.sg;
	tg = ci.tg;
	type = ci.type;
	valid = ci.valid;
	make_name();
}

cdup_item::cdup_item()
{
	type = CDUPLICATION;
	valid = true;
	make_name();
}

int cdup_item::make_name()
{
	name = "";
	adjacency::sort_adjacencies(s);

	char buf[1024];
	
	if(s.size() == 0)
	{
		assert(sg.size() == 1);
		sprintf(buf, "%d", sg.at(0)->x);
		name.append("C(");
		name.append(buf);
		name.append(")");
		return 0;
	}

	name = "C(";
	for(int i = 0; i < s.size(); i++)
	{
		sprintf(buf, "%d", s.at(i).e1.label());
		name.append(buf);
		name.append(",");
	}

	sprintf(buf, "%d", 0 - s.at(s.size() - 1).e2.label());
	name.append(buf);
	name.append(")");

	return 0;
}

int cdup_item::copy(const cdup_item & di)
{
	valid = di.valid;
	type = di.type;
	s = di.s;
	t = di.t;
	c = di.c;
	sg = di.sg;
	tg = di.tg;
	make_name();
	return 0;
}

/*
cdup_item::cdup_item(const PG & _y, const PG & _z, bool b)
{
	type = DUPLICATION;
	valid = true;

	if(b == true)
	{
		t.clear();
		tg.clear();
		gene * q = _y.first;
		while(true)
		{
			tg.push_back(q);
			if(q == _y.second) break;
			t.push_back(adjacency(PG(q, q->b)));
			q = q->b;
		}

		a1 = adjacency(PG(_z.first, _y.first));
		a2 = adjacency(PG(_y.second, _z.second));

		y = adjacency(_y);
		z = adjacency(_z);
	}
	else
	{
		t.clear();
		gene * q = _y.first;

		while(true)
		{
			tg.push_back(q);
			if(q == _y.second) break;
			t.push_back(adjacency(PG(q->a, q)));
			q = q->a;
		}

		a1 = adjacency(PG(_y.first, _z.first));
		a2 = adjacency(PG(_z.second, _y.second));

		y = adjacency(PG(_y.second, _y.first));
		z = adjacency(PG(_z.second, _z.first));
	}
}

int cdup_item::set_x(const PG & _x, bool b)
{
	vector<int> v;

	if(b == true)
	{
		s.clear();
		sg.clear();
		gene * q = _x.first;
		while(true)
		{
			sg.push_back(q);
			if(q == _x.second) break;
			s.push_back(adjacency(PG(q, q->b)));
			q = q->b;
		}
		x = adjacency(_x);
	}
	else
	{
		sg.clear();
		s.clear();
		gene * q = _x.first;
		while(true)
		{
			sg.push_back(q);
			if(q == _x.second) break;
			s.push_back(adjacency(PG(q->a, q)));
			q = q->a;
		}
		x = adjacency(PG(_x.second, _x.first));
	}

	assert(s.size() == t.size());

	make_name();

	return 0;
}

int cdup_item::set_z(const PG & _z, bool b)
{
	if(b == true) z = adjacency(_z);
	else z = adjacency(PG(_z.second, _z.first));
	return 0;
}
*/

int cdup_item::reverse_operate(genome_base & gm) const
{
	return 0;
}

int cdup_item::operate(genome_base & gm) const
{
	return 0;
}

int cdup_item::find_gene(gene * g)
{
	for(int i = 0; i < sg.size(); i++)
	{
		if(g == sg.at(i)) return i;
	}
	return -1;
}

int cdup_item::find_adjacency(const adjacency & a)
{
	for(int i = 0; i < s.size(); i++)
	{
		if(a == s.at(i)) return i + 1;
	}
	return -1;
}

bool cdup_item::overlap(const cdup_item & di)
{
	for(int i = 0; i < di.sg.size(); i++)
	{
		for(int j = 0; j < sg.size(); j++)
		{
			if(di.sg.at(i) == sg.at(j)) return true;
		}
	}
	return false;
}

int cdup_item::print() const
{
	printf("CDUPLICATION  : ");
	return 0;
}

/*
int cdup_item::rigid_compare(const cdup_item & x) const
{
	assert(x.sg.size() == x.tg.size());
	assert(sg.size() == tg.size());
	
	bool f = true;
	for(int i = 0; i < x.tg.size(); i++)
	{
		bool b = false;
		for(int k = 0; k < tg.size(); k++)
		{
			if(x.tg.at(i) == tg.at(k) && x.sg.at(i) == sg.at(k))
			{
				b = true;
				break;
			}
		}
		if(b == false)
		{
			f = false;
			break;
		}
	}

	if(f == true)
	{
		if(x.sg.size() == sg.size()) return EQUAL;
		else return SUBSET;
	}

	for(int i = 0; i < x.tg.size(); i++)
	{
		bool b = false;
		for(int k = 0; k < sg.size(); k++)
		{
			if(x.tg.at(i) == sg.at(k) && x.sg.at(i) == tg.at(k))
			{
				b = true;
				break;
			}
		}
		if(b == false) return OTHER;
	}

	if(x.tg.size() == sg.size()) return REVERSE_EQUAL;
	else return REVERSE_SUBSET;
}

int cdup_item::compare(const cdup_item & x) const
{
	assert(x.sg.size() == x.tg.size());
	assert(sg.size() == tg.size());

	
	bool f = true;
	for(int i = 0; i < x.tg.size(); i++)
	{
		bool b = false;
		for(int k = 0; k < tg.size(); k++)
		{
			if(x.tg.at(i) == tg.at(k))
			{
				b = true;
				break;
			}
		}
		if(b == false)
		{
			f = false;
			break;
		}
	}

	if(f == true)
	{
		if(x.sg.size() == sg.size()) return EQUAL;
		else return SUBSET;
	}

	for(int i = 0; i < x.tg.size(); i++)
	{
		bool b = false;
		for(int k = 0; k < sg.size(); k++)
		{
			if(x.tg.at(i) == sg.at(k))
			{
				b = true;
				break;
			}
		}
		if(b == false) return OTHER;
	}

	if(x.tg.size() == sg.size()) return REVERSE_EQUAL;
	else return REVERSE_SUBSET;
}
*/

int cdup_item::check_in_active(gene * g)
{
	for(int i = 0; i < sg.size(); i++)
	{
		if(g == sg.at(i)) return NONACTIVE_SOURCE;
	}
	for(int i = 0; i < tg.size(); i++)
	{
		if(g == tg.at(i)) return NONACTIVE_TARGET;
	}
	return NONACTIVE;
}

int cdup_item::check_out_active(gene * g)
{
	for(int i = 0; i < sg.size(); i++)
	{
		if(g == sg.at(i)) return NONACTIVE_SOURCE;
	}
	return NONACTIVE;
}

int cdup_item::check_in_active(const adjacency & a)
{
	for(int i = 0; i < s.size(); i++)
	{
		if(a == s.at(i)) return NONACTIVE_SOURCE;
	}
	for(int i = 0; i < t.size(); i++)
	{
		if(a == t.at(i)) return NONACTIVE_TARGET;
	}

	if(a == c) return NONACTIVE_TARGET;
	return NONACTIVE;
}

int cdup_item::check_out_active(const adjacency & a)
{
	for(int i = 0; i < s.size(); i++)
	{
		if(a == s.at(i)) return NONACTIVE_SOURCE;
	}
	return NONACTIVE;
}

pair<item*, item*> dup_item::split() const
{
	adjacency c = a1.diff(a2, z);
	item * di = new dcj_item(a1, a2, z, c);
	item * ci = new cdup_item(s, t, c, sg, tg);
	return pair<item*, item*>(di, ci);
}

dup_item::dup_item()
{
	type = DUPLICATION;
	valid = true;
}

int dup_item::make_name()
{
	name = "";
	adjacency::sort_adjacencies(s);

	char buf[1024];
	
	if(s.size() == 0)
	{
		assert(sg.size() == 1);
		sprintf(buf, "%d", sg.at(0)->x);
		name.append("D(");
		name.append(buf);
		name.append(")");
		return 0;
	}

	name = "D(";
	for(int i = 0; i < s.size(); i++)
	{
		sprintf(buf, "%d", s.at(i).e1.label());
		name.append(buf);
		name.append(",");
	}

	sprintf(buf, "%d", 0 - s.at(s.size() - 1).e2.label());
	name.append(buf);
	name.append(")");

	return 0;
}

int dup_item::copy(const dup_item & di)
{
	valid = di.valid;
	s = di.s;
	t = di.t;
	a1 = di.a1;
	a2 = di.a2;
	z = di.z;
	y = di.y;
	sg = di.sg;
	tg = di.tg;
	make_name();
	return 0;
}

dup_item::dup_item(const PG & _y, const PG & _z, bool b)
{
	type = DUPLICATION;
	valid = true;

	if(b == true)
	{
		t.clear();
		tg.clear();
		gene * q = _y.first;
		while(true)
		{
			tg.push_back(q);
			if(q == _y.second) break;
			t.push_back(adjacency(PG(q, q->b)));
			q = q->b;
		}

		a1 = adjacency(PG(_z.first, _y.first));
		a2 = adjacency(PG(_y.second, _z.second));

		y = adjacency(_y);
		z = adjacency(_z);
	}
	else
	{
		t.clear();
		gene * q = _y.first;

		while(true)
		{
			tg.push_back(q);
			if(q == _y.second) break;
			t.push_back(adjacency(PG(q->a, q)));
			q = q->a;
		}

		a1 = adjacency(PG(_y.first, _z.first));
		a2 = adjacency(PG(_z.second, _y.second));

		y = adjacency(PG(_y.second, _y.first));
		z = adjacency(PG(_z.second, _z.first));
	}
}


int dup_item::set_x(const PG & _x, bool b)
{
	vector<int> v;

	if(b == true)
	{
		s.clear();
		sg.clear();
		gene * q = _x.first;
		while(true)
		{
			sg.push_back(q);
			if(q == _x.second) break;
			s.push_back(adjacency(PG(q, q->b)));
			q = q->b;
		}
		//x = adjacency(_x);
	}
	else
	{
		sg.clear();
		s.clear();
		gene * q = _x.first;
		while(true)
		{
			sg.push_back(q);
			if(q == _x.second) break;
			s.push_back(adjacency(PG(q->a, q)));
			q = q->a;
		}
		//x = adjacency(PG(_x.second, _x.first));
	}

	assert(s.size() == t.size());

	make_name();

	return 0;
}

int dup_item::set_z(const PG & _z, bool b)
{
	if(b == true) z = adjacency(_z);
	else z = adjacency(PG(_z.second, _z.first));
	return 0;
}

int dup_item::reverse_operate(genome_base & gm) const
{
	assert(z.available() == true);

	for(int i = 0; i < s.size(); i++)
	{
		assert(s.at(i).available() == true);
	}
	
	adjacency zz = z;

	adjacency yy = a1.diff(a2, zz);
	yy.e1.b = yy.e1.b ? false : true;
	yy.e2.b = yy.e2.b ? false : true;

	if(zz.direction() == false) zz.exchange();
	if(yy.direction() == false) yy.exchange();

	adjacency o(zz.e1, yy.e1);
	o.e2.b = o.e2.b ? false : true;

	if(o != a1 && o != a2)
	{
		chrm::reverse(PG(yy.e1.g, yy.e2.g));
		yy.exchange();
	}

	yy.e1.g->a = NULL;
	yy.e2.g->b = NULL;
	operation * op = new insertion(zz.e2.g, yy.e1.g, yy.e2.g);

	gm.operate(op);
	delete op;

	for(int i = 0; i < t.size(); i++)
	{
		assert(t.at(i).available() == true);
	}

	for(int i = 0; i < s.size(); i++)
	{
		assert(s.at(i).available() == true);
	}

	assert(a1.available() == true);
	assert(a2.available() == true);

	return 0;
}

int dup_item::operate(genome_base & gm) const
{
	for(int i = 0; i < s.size(); i++)
	{
		assert(s.at(i).available() == true);
	}

	for(int i = 0; i < t.size(); i++)
	{
		assert(t.at(i).available() == true);
	}

	assert(a1.available() == true);
	assert(a2.available() == true);

	adjacency zz = z;
	if(zz.direction() == false) zz.exchange();

	operation * op = new deletion(zz.e1.g->b, zz.e2.g->a, false);
	gm.operate(op);
	delete op;

	assert(z.available() == true);
	return 0;
}

int dup_item::find_gene(gene * g)
{
	for(int i = 0; i < sg.size(); i++)
	{
		if(g == sg.at(i)) return i;
	}

	return -1;
}

int dup_item::find_adjacency(const adjacency & a)
{
	if(a == z) return 0;

	for(int i = 0; i < s.size(); i++)
	{
		if(a == s.at(i)) return i + 1;
	}

	return -1;
}

/*
bool dup_item::find_target_gene(gene * g)
{
	for(int j = 0; j < tg.size(); j++)
	{
		if(g == tg.at(j)) return true;
		if(g == tg.at(j)) return true;
	}
	return false;
}

bool dup_item::depend(const dup_item & di)
{
	for(int i = 0; i < di.sg.size(); i++)
	{
		for(int j = 0; j < tg.size(); j++)
		{
			if(di.sg.at(i) == tg.at(j)) return true;
		}
	}
	return false;
}
*/

bool dup_item::overlap(const dup_item & di)
{
	for(int i = 0; i < di.sg.size(); i++)
	{
		for(int j = 0; j < sg.size(); j++)
		{
			if(di.sg.at(i) == sg.at(j)) return true;
		}
	}
	return false;
}

int dup_item::print() const
{
	for(int i = 0; i < s.size(); i++)
	{
		printf("source %3d: ", i);
		s.at(i).print();
	}
	for(int i = 0; i < t.size(); i++)
	{
		printf("target %3d: ", i);
		t.at(i).print();
	}
	printf("context   : ");
	z.print();

	return 0;
}

int dup_item::rigid_compare(const dup_item & x) const
{
	assert(x.sg.size() == x.tg.size());
	assert(sg.size() == tg.size());
	
	bool f = true;
	for(int i = 0; i < x.tg.size(); i++)
	{
		bool b = false;
		for(int k = 0; k < tg.size(); k++)
		{
			if(x.tg.at(i) == tg.at(k) && x.sg.at(i) == sg.at(k))
			{
				b = true;
				break;
			}
		}
		if(b == false)
		{
			f = false;
			break;
		}
	}

	if(f == true)
	{
		if(x.sg.size() == sg.size()) return EQUAL;
		else return SUBSET;
	}

	for(int i = 0; i < x.tg.size(); i++)
	{
		bool b = false;
		for(int k = 0; k < sg.size(); k++)
		{
			if(x.tg.at(i) == sg.at(k) && x.sg.at(i) == tg.at(k))
			{
				b = true;
				break;
			}
		}
		if(b == false) return OTHER;
	}

	if(x.tg.size() == sg.size()) return REVERSE_EQUAL;
	else return REVERSE_SUBSET;
}

int dup_item::compare(const dup_item & x) const
{
	assert(x.sg.size() == x.tg.size());
	assert(sg.size() == tg.size());

	
	bool f = true;
	for(int i = 0; i < x.tg.size(); i++)
	{
		bool b = false;
		for(int k = 0; k < tg.size(); k++)
		{
			if(x.tg.at(i) == tg.at(k))
			{
				b = true;
				break;
			}
		}
		if(b == false)
		{
			f = false;
			break;
		}
	}

	if(f == true)
	{
		if(x.sg.size() == sg.size()) return EQUAL;
		else return SUBSET;
	}

	for(int i = 0; i < x.tg.size(); i++)
	{
		bool b = false;
		for(int k = 0; k < sg.size(); k++)
		{
			if(x.tg.at(i) == sg.at(k))
			{
				b = true;
				break;
			}
		}
		if(b == false) return OTHER;
	}

	if(x.tg.size() == sg.size()) return REVERSE_EQUAL;
	else return REVERSE_SUBSET;
}

int dup_item::check_in_active(gene * g)
{
	for(int i = 0; i < sg.size(); i++)
	{
		if(g == sg.at(i)) return NONACTIVE_SOURCE;
	}
	for(int i = 0; i < tg.size(); i++)
	{
		if(g == tg.at(i)) return NONACTIVE_TARGET;
	}
	return NONACTIVE;
}

int dup_item::check_out_active(gene * g)
{
	for(int i = 0; i < sg.size(); i++)
	{
		if(g == sg.at(i)) return NONACTIVE_SOURCE;
	}
	return NONACTIVE;
}

int dup_item::check_in_active(const adjacency & a)
{
	for(int i = 0; i < s.size(); i++)
	{
		if(a == s.at(i)) return NONACTIVE_SOURCE;
	}
	for(int i = 0; i < t.size(); i++)
	{
		if(a == t.at(i)) return NONACTIVE_TARGET;
	}

	if(a == a1 || a == a2) return ACTIVE;
	return NONACTIVE;
}

int dup_item::check_out_active(const adjacency & a)
{
	for(int i = 0; i < s.size(); i++)
	{
		if(a == s.at(i)) return NONACTIVE_SOURCE;
	}

	if(a == z) return ACTIVE;
	return NONACTIVE;
}

dcj_item::dcj_item(const dcj_item & di)
{
	a1 = di.a1;
	a2 = di.a2;
	o1 = di.o1;
	o2 = di.o2;
	type = di.type;
	valid = di.valid;
	make_name();
}

int dcj_item::copy(const dcj_item & di)
{
	type = di.type;
	valid = di.valid;
	a1 = di.a1;
	a2 = di.a2;
	o1 = di.o1;
	o2 = di.o2;
	make_name();
	return 0;
}

dcj_item::dcj_item(const adjacency & _a1, const adjacency & _a2, const adjacency & _o1, const adjacency & _o2) 
{
	a1 = _a1;
	a2 = _a2;
	o1 = _o1;
	o2 = _o2;
	type = DCJ;
	valid = true;
	make_name();
}

dcj_item::dcj_item(const adjacency & _a1, const adjacency & _a2, bool d)
	: a1(_a1), a2(_a2)
{
	type = DCJ;
	valid = true;

	if(d == true)
	{
		o1 = adjacency(a1.e1, a2.e2);
		o2 = adjacency(a2.e1, a1.e2);
	}
	else
	{
		o1 = adjacency(a1.e1, a2.e1);
		o2 = adjacency(a1.e2, a2.e2);
	}

	make_name();
}

int dcj_item::make_name()
{
	name = "J";
	return 0;
}

dcj_item::dcj_item(dcj * opx)
{
	type = DCJ;
	valid = true;
	//d = opx->d;

	a1 = adjacency(opx->pos1);
	a2 = adjacency(opx->pos2);

	if(opx->d == true)
	{
		o1 = adjacency(a1.e1, a2.e2);
		o2 = adjacency(a1.e2, a2.e1);
	}
	else
	{
		o1 = adjacency(a1.e1, a2.e1);
		o2 = adjacency(a1.e2, a2.e2);
	}

	make_name();
}

int dcj_item::reverse_operate(genome_base & gm) const
{
	assert(o1.available() == true);
	assert(o2.available() == true);

	adjacency oo1 = o1;
	adjacency oo2 = o2;

	if(oo1.direction() == false) oo1.exchange();
	if(oo2.direction() == false) oo2.exchange();

	assert(oo1.direction() == true);
	assert(oo2.direction() == true);

	PG pg1 = PG(oo1.e1.g, oo1.e2.g);
	PG pg2 = PG(oo2.e1.g, oo2.e2.g);

	assert(pg1.first->b == pg1.second);
	assert(pg2.first->b == pg2.second);

	operation * opx;

	if(a1 == adjacency(oo1.e1, oo2.e1) || a1 == adjacency(oo1.e2, oo2.e2))
	{
		opx = new dcj(pg1, pg2, false);
	}
	else if(a1 == adjacency(oo1.e1, oo2.e2) || a1 == adjacency(oo1.e2, oo2.e1))
	{
		opx = new dcj(pg1, pg2, true);
	}
	else
	{
		assert(1 == 0);
	}

	gm.operate(opx);
	delete opx;

	assert(a1.available() == true);
	assert(a2.available() == true);

	return 0;
}

int dcj_item::operate(genome_base & gm) const
{
	assert(a1.available() == true);
	assert(a2.available() == true);

	adjacency aa1 = a1;
	adjacency aa2 = a2;

	if(aa1.direction() == false) aa1.exchange();
	if(aa2.direction() == false) aa2.exchange();

	assert(aa1.direction() == true);
	assert(aa2.direction() == true);

	PG pg1 = PG(aa1.e1.g, aa1.e2.g);
	PG pg2 = PG(aa2.e1.g, aa2.e2.g);

	assert(pg1.first->b == pg1.second);
	assert(pg2.first->b == pg2.second);

	operation * opx;

	if(o1 == adjacency(aa1.e1, aa2.e1) || o1 == adjacency(aa1.e2, aa2.e2))
	{
		opx = new dcj(pg1, pg2, false);
	}
	else if(o1 == adjacency(aa1.e1, aa2.e2) || o1 == adjacency(aa1.e2, aa2.e1))
	{
		opx = new dcj(pg1, pg2, true);
	}
	else
	{
		assert(1 == 0);
	}

	gm.operate(opx);
	delete opx;

	assert(o1.available() == true);
	assert(o2.available() == true);

	return 0;
}

int dcj_item::print() const
{
	printf("DCJ          : ");	
	printf("%s + %s -> %s + %s\n", a1.label_str().c_str(), a2.label_str().c_str(),
			o1.label_str().c_str(), o2.label_str().c_str());
	return 0;
}

int dcj_item::find_adjacency(const adjacency & a)
{
	if(a == o1) return 0;
	if(a == o2) return 1;
	return -1;
}

int dcj_item::find_gene(gene * g)
{
	return -1;
}

int dcj_item::reverse()
{
	adjacency x1 = a1;
	adjacency x2 = a2;
	a1 = o1;
	a2 = o2;
	o1 = x1;
	o2 = x2;
	return 0;
}

/*
dcj_item dcj_item::reverse()
{
	dcj_item di(*this);
	di.reverse();
	return di;
}
*/

int dcj_item::compare(const dcj_item & x) const
{
	assert(a1 != a2);
	assert(o1 != o2);
	assert(x.a1 != x.a2);
	assert(x.o1 != x.o2);

	bool equal = true;
	if(x.a1 != a1 && x.a1 != a2) equal = false;
	if(x.a2 != a1 && x.a2 != a2) equal = false;
	if(x.o1 != o1 && x.o1 != o2) equal = false;
	if(x.o2 != o1 && x.o2 != o2) equal = false;
	if(equal == true) return EQUAL;

	bool reverse_equal = true;
	if(x.a1 != o1 && x.a1 != o2) reverse_equal = false;
	if(x.a2 != o1 && x.a2 != o2) reverse_equal = false;
	if(x.o1 != a1 && x.o1 != a2) reverse_equal = false;
	if(x.o2 != a1 && x.o2 != a2) reverse_equal = false;
	if(reverse_equal == true) return REVERSE_EQUAL;

	return OTHER;
}

int dcj_item::check_in_active(gene * g)
{
	return NONACTIVE;
}

int dcj_item::check_out_active(gene * g)
{
	return NONACTIVE;
}

int dcj_item::check_in_active(const adjacency & a)
{
	if(a == a1 || a == a2) return ACTIVE;
	return NONACTIVE;
}

int dcj_item::check_out_active(const adjacency & a)
{
	if(a == o1 || a == o2) return ACTIVE;
	return NONACTIVE;
}
