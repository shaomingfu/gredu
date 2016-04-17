#include "operation.h"
#include "circular_chrm.h"
#include "common.h"

#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <iostream>

// both start and end position are inclusive

circular_chrm::circular_chrm()
	: chrm(CIRCULAR)
{
	p = NULL;
}

circular_chrm::circular_chrm(const chrm & lc)
	: chrm(lc)
{
	if(lc.p == NULL) return;
	assert(lc.p->a != NULL);
	assert(lc.p->b != NULL);

	//p = new gene(*(lc.p));
	p = new gene(lc.p->x, this);

	assert(p->a == NULL);
	assert(p->b == NULL);

	gene * x = lc.p;
	gene * y = p;

	while(x->b != lc.p)
	{
		//y->b = new gene(*(x->b));
		y->b = new gene(x->b->x, this);
		y->b->a = y;

		x = x->b;
		y = y->b;
	}

	assert(x->b == lc.p);
	y->b = p;
	p->a = y;
}

circular_chrm::~circular_chrm()
{
	destroy();
}

int circular_chrm::destroy()
{
	if(p == NULL) return 0;

	assert(p->a != NULL);
	assert(p->b != NULL);

	gene * q = p;
	while(q->b != p)
	{
		gene * x = q->b;
		delete q;
		q = x;
	}

	assert(q->b == p);
	delete q;

	p = NULL;

	return 0;
}

int circular_chrm::size()
{
	if(p == NULL) return 0;

	gene * q = p;

	int x = 0;
	while(true)
	{
		assert(p->x != 0);
		x++;
		if(q->b == p) break;
		q = q->b;
	}

	return x;
}

gene * circular_chrm::back() const
{
	if(p == NULL) return NULL;

	assert(p->a != NULL);
	assert(p->b != NULL);

	return p->a;
}

int circular_chrm::regular_print(const PG & pg)
{
	gene * q = pg.first;
	vector<gene*> v;
	while(true)
	{
		v.push_back(q);
		if(q == pg.second) break;
		q = q->b;
	}

	regular_print(v);

	return 0;
}

int circular_chrm::regular_print(gene * g)
{
	if(p == NULL) 
	{
		printf("[circular] [empty]\n");
		return 0;
	}
	
	printf("[circular]    ");

	gene * q = p;
	assert(q->b != NULL);

	while(q->b != p)
	{
		assert(q->ch == this);

		if(q == g) printf("%3d ", g->x);
		else printf("    ");
		q = q->b;
	}

	if(q == g) printf("%3d ", g->x);
	else printf("    ");

	return 0;
}


int circular_chrm::regular_print(const PG & pg, const PI & pi)
{
	if(p == NULL) 
	{
		printf("[circular] [empty]\n");
		return 0;
	}
	
	printf("[circular]    ");

	gene * q = p;
	assert(q->b != NULL);

	while(q->b != p)
	{
		assert(q->ch == this);

		if(q == pg.first) printf("%3d ", pi.first);
		else if(q == pg.second) printf("%3d ", pi.second);
		else printf("    ");

		q = q->b;
	}

	if(q == pg.first) printf("%3d ", pi.first);
	else if(q == pg.second) printf("%3d ", pi.second);
	else printf("    ");

	return 0;
}

int circular_chrm::regular_print(const vector<gene*> & list)
{
	if(p == NULL) 
	{
		printf("[circular] [empty]\n");
		return 0;
	}
	
	printf("[circular]    ");

	gene * q = p;
	assert(q->b != NULL);

	while(q->b != p)
	{
		assert(q->ch == this);

		bool b = false;
		for(int i = 0; i < list.size(); i++)
		{
			if(list.at(i) == q) b = true; 
		}

		if(b == true) printf("%3d ", q->x);
		else printf("    ");

		q = q->b;
	}

	bool b = false;
	for(int i = 0; i < list.size(); i++)
	{
		if(list.at(i) == q) b = true; 
	}

	if(b == true) printf("%3d \n", q->x);
	else printf("    \n");

	return 0;
}

int circular_chrm::regular_print()
{
	if(p == NULL) 
	{
		printf("[circular] [empty]\n");
		return 0;
	}
	
	printf("[circular]    ");

	gene * q = p;
	assert(q->b != NULL);

	while(q->b != p)
	{
		assert(q->ch == this);
		printf("%3d ", q->x);
		q = q->b;
	}
	printf("%3d \n", q->x);

	return 0;
}

int circular_chrm::blank_print()
{
	if(p == NULL) 
	{
		return 0;
	}
	
	gene * q = p;
	assert(q->b != NULL);

	while(q->b != p)
	{
		assert(q->ch == this);
		cout<<q->x<<" ";
		q = q->b;
	}
	cout<<q->x;
	cout<<endl;

	return 0;
}


int circular_chrm::print()
{
	if(p == NULL) 
	{
		cout<<"[circular] [empty]\n";
		return 0;
	}
	
	cout<<"[circular]    ";

	gene * q = p;
	assert(q->b != NULL);

	while(q->b != p)
	{
		assert(q->ch == this);
		cout<<q->x<<" ";
		q = q->b;
	}
	cout<<q->x;
	cout<<endl;

	return 0;
}

int circular_chrm::print_address()
{
	if(p == NULL) 
	{
		cout<<"[circular] [empty]\n";
		return 0;
	}
	
	cout<<"[circular]("<<this<<", "<<p<<"):  ";

	gene * q = p;
	assert(q->b != NULL);

	while(q->b != p)
	{
		assert(q->ch == this);
		cout<<q->x<<"("<<q->ch<<")  ";
		q = q->b;
	}
	cout<<q->x<<"("<<q->ch<<")";
	cout<<endl;

	return 0;
}

int circular_chrm::reverse_print()
{
	gene * b = back();

	if(b == NULL) 
	{
		cout<<"[circular:r] [empty]\n";
		return 0;
	}
	
	assert(b->b != NULL);

	cout<<"[circular:r]  ";

	gene * q = b;
	while(q->a != b)
	{
		assert(q->ch == this);
		cout<<0 - q->x<<" ";
		q = q->a;
	}
	cout<<0 - q->x;
	cout<<endl;

	return 0;
}

int circular_chrm::insert_front(int x)
{
	gene * g = new gene(x);

	if(p == NULL)
	{
		p = g;
		p->a = p;
		p->b = p;
	}
	else
	{
		g->b = p;
		g->a = p->a;
		p->a->b = g;
		p->a = g;
		p = g;
	}

	return 0;
}

int circular_chrm::init(const vector<int> & s)
{
	destroy();

	assert(p == NULL);

	for(int i = 0; i < s.size(); i++)
	{
		insert_front(s.at(i));
	}

	return 0;
}

int circular_chrm::random_init(int len, int alphabet_size)
{
	destroy();

	assert(p == NULL);

	for(int i = 0; i < len; i++)
	{
		insert_front( rand() % alphabet_size + 1 );
	}

	return 0;
}

int circular_chrm::reverse()
{
	if(p == NULL) return 0;

	gene * q = p;
	gene * b = back();

	while(true)
	{
		q->reverse();
		if(q == b) break;
		q = q->a;
	}

	return 0;
}


int circular_chrm::enumerate_reverse_duplication(chrm * ch, vector<operation*> & v, int max_length)
{
	assert(false);
	return 0;
}

int circular_chrm::enumerate_reverse_duplication(vector<operation*> & v, int max_length)
{
	int n = size();

	if(max_length > n) max_length = n;
	for(int l = 1; l <= max_length; l++)
	{
		gene * pi = p;
		for(int i = 0; i < n - l; i++)
		{
			gene * pj = chrm::advance(pi, l);
			for(int j = i + l; j < n; j++)
			{
				// check
				bool b = true;
				gene * pik = pi;
				gene * pjk = pj;
				for(int k = 0; k < l; k++)
				{
					if(pik->x != pjk->x)
					{
						b = false;
						break;
					}
					pik = pik->b;
					pjk = pjk->b;
				}
				
				if(b==true)
				{
					v.push_back(new duplication(pi, pj, advance(pj, l - 1)));
					v.push_back(new duplication(pj, pi, advance(pi, l - 1)));
				}
				pj = pj->b;
			}
			pi = pi->b;
		}
	}
	return 0;
}


int circular_chrm::enumerate_reverse_iduplication(chrm * ch, vector<operation*> & v, int max_length)
{
	assert(false);
	return 0;
}

int circular_chrm::enumerate_reverse_iduplication(vector<operation*> & v, int max_length)
{
	int n = size();

	if(max_length > n) max_length = n;

	for(int l = 1; l <= max_length; l++)
	{
		gene * pi = p;
		for(int i = 0; i < n - l; i++)
		{
			gene * pj = chrm::advance(pi, l);
			for(int j = i + l; j < n; j++)
			{
				// check
				bool b = true;
				gene * pik = pi;
				gene * pjk = chrm::advance(pj, l - 1);
				for(int k = 0; k < l; k++)
				{
					if(pik->x != 0 - pjk->x)
					{
						b = false;
						break;
					}
					pik = pik->b;
					pjk = pjk->a;
				}
				
				if(b==true)
				{
					v.push_back(new iduplication(pi, pj, advance(pj, l - 1)));
					v.push_back(new iduplication(pj, pi, advance(pi, l - 1)));
				}
				pj = pj->b;
			}
			pi = pi->b;
		}
	}
	return 0;
}


int circular_chrm::count_duplication(int max_length)
{
	int n = size();
	if(max_length > n) max_length = n;
	return n * n * max_length;
}

gene * circular_chrm::at(int index)
{
	assert(index >= 0 && index < size());
	return chrm::advance(p, index);
}

int circular_chrm::index(gene * g)
{
	int x = 0;
	gene * q = p;
	while(true)
	{
		if(g == q) return x;
		assert(q->b != p);
		x++;
		q = q->b;
		if(q == p) return -1;
	}
	//assert(1 == 0);
	return -1;
}

bool circular_chrm::is_equal(const chrm * ch)
{
	assert(1 == 0);
	return false;
}

