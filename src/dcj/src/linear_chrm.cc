#include "operation.h"
#include "linear_chrm.h"
#include "common.h"


#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <iostream>


linear_chrm::linear_chrm()
	: chrm(LINEAR)
{
	p = new gene(0, this);
	p->b = new gene(0, this);
	p->b->a = p;
}

/*
linear_chrm::linear_chrm(gene * end1, gene * end2)
	: chrm(LINEAR)
{
	assert(end1 != NULL);
	assert(end2 != NULL);
	assert(end1->x == 0);
	assert(end2->x == 0);

	p = end1;
	p->b = end2;
	p->a = NULL;
	p->b->b = NULL;
	p->ch = this;
	p->b->ch = this;
}
*/

linear_chrm::linear_chrm(const chrm & lc)
	: chrm(lc)
{
	assert(lc.p != NULL);
	assert(lc.p->a == NULL);
	assert(lc.p->x == 0);

	//p = new gene(*(lc.p));
	p = new gene(lc.p->x, this);
	assert(p->a == NULL);

	gene * x = lc.p;
	gene * y = p;

	while(x->b != NULL)
	{
		//y->b = new gene(*(x->b));
		y->b = new gene(x->b->x, this);
		y->b->a = y;

		x = x->b;
		y = y->b;
	}

	assert(x->b == NULL);
	assert(y->b == NULL);
}

linear_chrm::~linear_chrm()
{
	destroy();
}

int linear_chrm::destroy()
{
	assert(p != NULL);
	assert(p->a == NULL);

	gene * q = p;
	while(q != NULL)
	{
		gene * x = q->b;
		delete q;
		q = x;
	}

	p = NULL;
	return 0;
}

int linear_chrm::size()
{
	assert(p != NULL);

	int x = 0;
	gene * q = p;
	while(true)
	{
		if(q->x != 0) x++;
		if(q->b == NULL)
		{
			assert(q->x == 0);
			break;
		}
		q = q->b;
	}

	return x;
}

gene * linear_chrm::back() const
{
	assert(p != NULL);

	gene * q = p;

	while(q->b != NULL) q = q->b;

	return q;
}

int linear_chrm::print()
{
	assert(p != NULL);

	cout<<"[linear]      ";
	gene * q = p;
	assert(q->a == NULL);

	while(q != NULL)
	{
		assert(q->ch == this);
		cout<<q->x<<" ";
		q = q->b;
	}
	cout<<endl;

	return 0;
}

int linear_chrm::blank_print()
{
	assert(p != NULL);

	gene * q = p;
	assert(q->a == NULL);

	while(q != NULL)
	{
		assert(q->ch == this);
		if(q->x != 0) cout<<q->x<<" ";
		q = q->b;
	}
	cout<<endl;

	return 0;
}

int linear_chrm::regular_print()
{
	assert(p != NULL);

	printf("[linear]      ");
	gene * q = p;
	assert(q->a == NULL);

	while(q != NULL)
	{
		assert(q->ch == this);
		printf("%3d ", q->x);
		q = q->b;
	}
	printf("\n");

	return 0;
}

int linear_chrm::regular_print(gene * g)
{
	printf("[linear]      ");

	gene * q = p;
	assert(q->b != NULL);

	while(q != NULL)
	{
		assert(q->ch == this);

		if(q == g) printf("%3d ", g->x);
		else printf("    ");

		q = q->b;
	}
	printf("\n");

	return 0;
}

int linear_chrm::regular_print(const PG & pg, const PI & pi)
{
	printf("[linear]      ");

	gene * q = p;
	assert(q->b != NULL);

	while(q != NULL)
	{
		assert(q->ch == this);

		if(q == pg.first) printf("%3d ", pi.first);
		else if(q == pg.second) printf("%3d ", pi.second);
		else printf("    ");

		q = q->b;
	}
	printf("\n");

	return 0;
}

int linear_chrm::regular_print(const vector<gene*> & list)
{
	printf("[linear]      ");

	gene * q = p;
	assert(q->b != NULL);

	while(q != NULL)
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
	printf("\n");

	return 0;
}

int linear_chrm::regular_print(const PG & pg)
{
	assert(p != NULL);

	printf("              ");

	gene * q = p;
	assert(q->a == NULL);

	bool start = false;
	bool end = false;
	while(q != NULL)
	{
		assert(q->ch == this);

		if(q == pg.first) start = true;

		if(start == true && end == false) printf("%3d ", q->x);
		else printf("    ");

		if(q == pg.second) end = true;

		q = q->b;
	}
	printf("\n");

	return 0;
}

int linear_chrm::print_address()
{
	assert(p != NULL);

	cout<<"[linear]("<<this<<", "<<p<<"):  ";
	gene * q = p;
	assert(q->a == NULL);

	while(q != NULL)
	{
		assert(q->ch == this);
		cout<<q->x<<"("<<q->ch<<")  ";
		q = q->b;
	}
	cout<<endl;

	return 0;
}

int linear_chrm::reverse_print()
{
	gene * b = back();

	assert(b != NULL);

	assert(b->b == NULL);

	cout<<"[linear:r]    ";

	gene * q = b;
	while(q != NULL)
	{
		assert(q->ch == this);
		cout<<0 - q->x<<" ";
		q = q->a;
	}
	cout<<endl;

	return 0;
}

int linear_chrm::insert_front(int x)
{
	gene * g = new gene(x, this);

	g->a = p;
	p->b = g;

	g->b = p->b;
	p->b->a = g;

	return 0;
}

int linear_chrm::init(const vector<int> & s)
{
	destroy();

	assert(p == NULL);

	for(int i = 0; i < s.size(); i++)
	{
		insert_front(s.at(i));
	}

	return 0;
}

int linear_chrm::random_init(int len, int alphabet_size)
{
	destroy();

	assert(p == NULL);

	for(int i = 0; i < len; i++)
	{
		insert_front( rand() % alphabet_size + 1 );
	}

	return 0;
}


int linear_chrm::reverse()
{
	gene * q = p;
	gene * b = back();
	while(true)
	{
		q->reverse();
		if(q == b) break;
		q = q->a;
	}

	p = b;
	return 0;
}

int linear_chrm::count_duplication(int max_length)
{
	int n = size();

	if(max_length > n) max_length = n;

	int x = max_length;

	return (n + 1) * n * x - (n + 1) * x * (x - 1) / 2;
}

int linear_chrm::enumerate_reverse_duplication(chrm * ch, vector<operation*> & v, int max_length)
{
	int nx = size();
	int ny = ch->size();
	if(max_length > nx) max_length = nx;
	if(max_length > ny) max_length = ny;

	linear_chrm * lch = dynamic_cast<linear_chrm*>(ch);

	for(int l = 1; l <= max_length; l++)
	{
		gene * pi = p->b;
		for(int i = 0; i <= nx - l; i++)
		{
			//gene * pj = chrm::advance(pi, l);
			gene * pj = lch->p->b;
			for(int j = 0; j <= ny - l; j++)
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


int linear_chrm::enumerate_reverse_duplication(vector<operation*> & v, int max_length)
{
	int n = size();

	if(max_length > n) max_length = n;

	for(int l = 1; l <= max_length; l++)
	{
		gene * pi = p->b;
		for(int i = 0; i <= n - l; i++)
		{
			gene * pj = chrm::advance(pi, l);
			for(int j = i + l; j <= n - l; j++)
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

int linear_chrm::enumerate_reverse_iduplication(chrm * ch, vector<operation*> & v, int max_length)
{
	int nx = size();
	int ny = ch->size();
	if(max_length > nx) max_length = nx;
	if(max_length > ny) max_length = ny;

	linear_chrm * lch = dynamic_cast<linear_chrm*>(ch);

	for(int l = 1; l <= max_length; l++)
	{
		gene * pi = p->b;
		for(int i = 0; i <= nx - l; i++)
		{
			gene * pj = lch->p->b;
			for(int j = 0; j <= ny - l; j++)
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

int linear_chrm::enumerate_reverse_iduplication(vector<operation*> & v, int max_length)
{
	int n = size();

	if(max_length > n) max_length = n;

	for(int l = 1; l <= max_length; l++)
	{
		gene * pi = p->b;
		for(int i = 0; i <= n - l; i++)
		{
			gene * pj = chrm::advance(pi, l);
			for(int j = i + l; j <= n - l; j++)
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

gene * linear_chrm::at(int index)
{
	if(index == -1) return p;
	assert(index >= 0 && index <= size());
	return chrm::advance(p->b, index);
}

int linear_chrm::index(gene * g)
{
	if(g == p) return -1;

	gene * q = p->b;
	int x = 0;
	while(true)
	{
		if(q == g) return x;
		q = q->b;
		x++;
		if(q->b == NULL) 
		{
			assert(q->x == 0);
			return -1;
		}
	}
	//assert(1 == 0);
	return -1;
}

bool linear_chrm::is_equal(const chrm * ch)
{
	bool b1 = chrm::is_equal(PG(p, back()), PG(ch->p, ch->back()));
	if(b1 == true) return true;

	bool b2 = chrm::is_reverse_equal(PG(p, back()), PG(ch->p, ch->back()));
	if(b2 == true) return true;

	return false;
}
