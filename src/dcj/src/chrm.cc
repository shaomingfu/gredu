#include "chrm.h"
#include "common.h"

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cassert>

chrm::chrm(int _type)
	: type(_type)
{
	p = NULL;
}

chrm::chrm(const chrm & c)
{
	type = c.type;
	p = NULL;
}

chrm::~chrm()
{}

PG chrm::make_gene_list(const vector<int> & s, chrm * ch)
{
	if(s.size() == 0) return PG(NULL, NULL);

	gene * p = new gene(s.at(0), ch);
	assert(p->a == NULL);
	assert(p->b == NULL);

	gene * q = p;
	for(int i = 1; i < s.size(); i++)
	{
		q->b = new gene(s.at(i), ch);
		q->b->a = q;
		q = q->b;
	}

	assert(p->a == NULL);
	assert(q->b == NULL);
	return PG(p, q);
}

PG chrm::make_gene_list(const vector<int> & s, const vector<string> & ss, chrm * ch)
{
	assert(s.size() == ss.size());
	if(s.size() == 0) return PG(NULL, NULL);

	gene * p = new gene(s.at(0), ch, ss.at(0));
	assert(p->a == NULL);
	assert(p->b == NULL);

	gene * q = p;
	for(int i = 1; i < s.size(); i++)
	{
		q->b = new gene(s.at(i), ch, ss.at(i));
		q->b->a = q;
		q = q->b;
	}

	assert(p->a == NULL);
	assert(q->b == NULL);
	return PG(p, q);
}

PG chrm::duplicate(gene * head, gene * tail)
{
	assert(head != NULL);
	assert(tail != NULL);

	gene * x = new gene(*head);

	gene * p = x;
	gene * q = head;
	while(true)
	{
		if(q == tail) break;
		p->b = new gene(*(q->b));
		p->b->a = p;
		p = p->b;
		q = q->b;
	}

	assert(x->x == head->x);
	assert(p->x == tail->x);

	x->a = NULL;
	p->b = NULL;

	return PG(x, p);
}

PG chrm::iduplicate(gene * head, gene * tail)
{
	PG pg = duplicate(head, tail);

	gene * q = pg.first;

	while(true)
	{
		q->reverse();
		if(q == pg.second) break;
		q = q->a;
	}

	return PG(pg.second, pg.first);
}

int chrm::count_inversion()
{
	int n = size();
	return n * (n + 1);
}

int chrm::count_loss()
{
	return size();
}

bool chrm::order(gene * x, gene * y)
{
	assert(x != NULL);
	assert(y != NULL);

	assert(x->ch == y->ch);

	gene * q = x;

	while(true)
	{
		if(q == y) return true;
		q = q->b;
		if(q == NULL) break;
	}

	q = y;

	while(true)
	{
		if(q == x) return false;
		q = q->b;
		assert(q != NULL);
	}

	assert(1 == 0);

	return false;
}

gene * chrm::advance(gene * g, int steps)
{
	gene * x = g;
	for(int i = 0; i < steps; i++)
	{
		x = x->b;
	}

	return x;
}

int chrm::print(gene * head, gene * tail)
{
	gene * g = head;
	printf("(");
	while(true)
	{
		assert(g != NULL);
		if(g == tail)
		{
			printf("%d)", g->x);
			break;
		}
		else
		{
			printf("%d ", g->x);
			g = g->b;
		}
	}
	return 0;
}

int chrm::print(gene * head, int len)
{
	if(len == 0) return 0;

	gene * g = head;
	printf("(");
	for(int i = 0; i < len; i++)
	{
		if(i < len - 1)
		{
			printf("%d ", g->x);
			g = g->b;
		}
		else
		{
			printf("%d)", g->x);
			return 0;
		}
	}
	return 0;
}

PG chrm::sort(const PG & pg)
{
	if(pg.second->b == pg.first) return PG(pg.second, pg.first);

	assert(pg.first->b == pg.second);

	return pg;
}

bool chrm::is_equal(const PG& pg1, const PG& pg2)
{
	gene * p = pg1.first;
	gene * q = pg2.first;

	while(true)
	{
		if(p->x != q->x) return false;
		if(p == pg1.second)
		{
			if(q != pg2.second) return false;
			break;
		}
		p = p->b;
		q = q->b;
	}
	return true;
}

bool chrm::is_equal(const PG & pg, const vector<int> & v)
{
	gene * q = pg.first;
	int k = 0;
	while(true)
	{
		if(q->x != v.at(k)) return false;
		if(q == pg.second)
		{
			if(k != v.size() - 1) return false;
			break;
		}
		k++;
		q = q->b;
	}
	return true;
}

bool chrm::is_reverse_equal(const PG & pg1, const PG & pg2)
{
	gene * p = pg1.first;
	gene * q = pg2.second;

	while(true)
	{
		if(p->x != 0 - q->x) return false;
		if(p == pg1.second)
		{
			if(q != pg2.first) return false;
			else return true;
		}

		if(q == pg2.first)
		{
			if(p != pg1.second) return false;
			else return true;
		}

		p = p->b;
		q = q->a;
	}

	assert(1 == 0);
}

bool chrm::is_reverse_equal(const PG & pg, const vector<int> & v)
{
	vector<int> vv = v;
	invert(vv);
	return is_equal(pg, vv);
}

int chrm::assert_equal(const PG& pg1, const PG& pg2)
{
	assert( is_equal(pg1, pg2) == true );
	return 0;
}


bool chrm::is_overlap(const PG & pg1, const PG & pg2)
{
	gene * p = pg1.first;

	while(true)
	{
		if(p == pg2.first) return true;
		if(p == pg2.second) return true;
		if(p == pg1.second) break;
		p = p->b;
	}

	p = pg2.first;

	while(true)
	{
		if(p == pg1.first) return true;
		if(p == pg1.second) return true;
		if(p == pg2.second) break;
		p = p->b;
	}
	
	return false;
}

int chrm::reverse(const PG & pg)
{
	gene * q = pg.first;
	gene * b = pg.second;
	while(true)
	{
		q->reverse();
		if(q == b) break;
		q = q->a;
	}
	return 0;
}
