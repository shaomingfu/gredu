#include "genome_base.h"
#include "chrm.h"
#include "linear_chrm.h"
#include "circular_chrm.h"

#include <vector>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <iostream>

using namespace std;

genome_base::genome_base()
{
	//alphabet_size = 0;
}

/*
genome_base::genome_base(int _alphabet_size)
{
	alphabet_size = _alphabet_size;
}
*/

genome_base& genome_base::operator=(const genome_base & g)
{
	//alphabet_size = g.alphabet_size;

	for(int i = 0; i < chrms.size(); i++)
	{
		if(chrms.at(i) != NULL)
		{
			delete chrms.at(i);
		}
	}

	chrms.clear();

	for(int i = 0; i < g.chrms.size(); i++)
	{
		if(g.chrms.at(i)->type==LINEAR)
		{
			chrms.push_back( new linear_chrm( *(g.chrms.at(i)) ) );
		}
		else
		{
			chrms.push_back( new circular_chrm( *(g.chrms.at(i)) ) );
		}
	}

	return (*this);
}

genome_base::genome_base(const genome_base & g)
{
	//alphabet_size = g.alphabet_size;
	chrms.clear();
	for(int i = 0; i < g.chrms.size(); i++)
	{
		if(g.chrms.at(i)->type==LINEAR)
		{
			chrms.push_back( new linear_chrm( *(g.chrms.at(i)) ) );
		}
		else
		{
			chrms.push_back( new circular_chrm( *(g.chrms.at(i)) ) );
		}
	}
}

genome_base::~genome_base()
{
	for(int i = 0; i < chrms.size(); i++)
	{
		if(chrms.at(i) != NULL)
		{
			delete chrms.at(i);
		}
	}
}

int genome_base::operate(const PO & po)
{
	if(po.first != NULL) operate(po.first);
	if(po.second != NULL) operate(po.second);
	return 0;
}

int genome_base::operate(operation * op)
{
	if(op->type == INVERSION) do_inversion(op);
	else if(op->type == INSERTION) do_insertion(op);
	else if(op->type == DELETION) do_deletion(op);
	else if(op->type == CDUPLICATION) do_cduplication(op);
	else if(op->type == DUPLICATION) do_duplication(op);
	else if(op->type == IDUPLICATION) do_iduplication(op);
	else if(op->type == LOSS) do_loss(op);
	else if(op->type == CIRCULARIZATION) do_circularization(op);
	else if(op->type == EXCHANGE) do_exchange1(op);
	else if(op->type == DCJ) do_dcj(op);
	else if(op->type == REVERSION) do_reversion(op);
	else if(op->type == LINEARIZATION) do_linearization(op);
	else assert(1 == 0);
	return 0;
}

PO genome_base::reverse(operation * op) const
{
	if(op->type == DCJ) return reverse_dcj(op);
	else if(op->type == REVERSION) return reverse_reversion(op);
	else if(op->type == LINEARIZATION) return reverse_linearization(op);
	else if(op->type == CIRCULARIZATION) return reverse_circularization(op);
	else assert(1 == 0);
}

PO genome_base::reverse_circularization(operation * op) const
{
	assert(op->type == CIRCULARIZATION);
	const circularization * opx = dynamic_cast<circularization*>(op);
	return PO(new linearization(opx->head, opx->tail), NULL);
}

PO genome_base::reverse_linearization(operation * op) const
{
	assert(op->type == LINEARIZATION);
	const linearization * opx = dynamic_cast<linearization*>(op);
	return PO(new circularization(opx->head, opx->tail), NULL);
}

PO genome_base::reverse_reversion(operation * op) const
{
	assert(op->type == REVERSION);
	const reversion * opx = dynamic_cast<reversion *>(op);
	return PO(new reversion(opx->g), NULL);
}

PO genome_base::reverse_dcj(operation * op) const
{
	assert(op->type == DCJ);

	const dcj * opx = dynamic_cast<dcj *>(op);

	PG npos1 = opx->pos1;
	PG npos2 = opx->pos2;

	assert(npos1.first != NULL);
	assert(npos1.second != NULL);
	assert(npos2.first != NULL);
	assert(npos2.second != NULL);

	assert(npos1.first->ch == npos1.second->ch);
	assert(npos2.first->ch == npos2.second->ch);
	assert(npos1.first->b == npos1.second);
	assert(npos2.first->b == npos2.second);

	PO po(NULL, NULL);

	/*
	if(npos2.first->x == 0 && npos2.second->x == 0)
	{
		assert(npos2.first->a == NULL);
		assert(npos2.second->b == NULL);
		linearization * lop = new linearization(npos1.first, npos1.second, npos2.first, npos2.second);
		po.first = new dcj(PG(npos2.first, npos1.first), PG(npos1.second, npos2.second), true);
	}
	*/
	if(npos1.first->ch != npos2.first->ch)
	{
		if(opx->d == false)
		{
			npos2 = PG(npos2.second, npos2.first);
			po.second = new reversion(npos2.first);
		}

		po.first = new dcj(PG(npos1.first, npos2.second), PG(npos2.first, npos1.second), true);
	}
	else
	{
		bool t = chrm::order(npos1.first, npos2.first);
		if(t == false)
		{
			PG t = npos1;
			npos1 = npos2;
			npos2 = t;
		}

		if(opx->d == true)
		{
			po.first = new dcj(PG(npos1.first, npos2.second), PG(npos2.first, npos1.second), true);
		}
		else
		{
			po.first = new dcj(PG(npos1.first, npos2.first), PG(npos1.second, npos2.second), false);
		}
	}

	assert(po.first->type == DCJ);

	return po;
}

int genome_base::add_linear_chrm(const vector<int> & s)
{
	chrm * chr = new linear_chrm();
	chrms.push_back(chr);

	if(s.size() == 0) return 0;

	PG pg = chrm::make_gene_list(s, chr);

	operation * op = new insertion(chr->p->b, pg.first, pg.second);

	operate(op);

	delete op;

	return 0;
}

int genome_base::add_circular_chrm(const vector<int> & s)
{
	chrm * chr = new circular_chrm();
	chrms.push_back(chr);

	if(s.size() == 0) return 0;

	PG pg = chrm::make_gene_list(s, chr);

	operation * op = new insertion(chr->p, pg.first, pg.second);

	operate(op);

	delete op;

	return 0;
}

int genome_base::add_linear_chrm(const vector<int> & s, const vector<string> & ss)
{
	if(s.size() == 0) return 0;

	chrm * chr = new linear_chrm();
	chrms.push_back(chr);

	PG pg = chrm::make_gene_list(s, ss, chr);

	operation * op = new insertion(chr->p->b, pg.first, pg.second);

	operate(op);

	delete op;

	return 0;
}

int genome_base::add_circular_chrm(const vector<int> & s, const vector<string> & ss)
{
	if(s.size() == 0) return 0;

	chrm * chr = new circular_chrm();
	chrms.push_back(chr);

	PG pg = chrm::make_gene_list(s, ss, chr);

	operation * op = new insertion(chr->p, pg.first, pg.second);

	operate(op);

	delete op;

	return 0;
}

int genome_base::do_exchange1(operation * op)
{
	assert(op->type == EXCHANGE);
	const exchange1 * opx = dynamic_cast<exchange1 *>(op);

	/*
	printf("=================\n");
	print();
	printf("-----------------\n");
	opx->pos1.first->ch->print();
	printf("-----------------\n");
	opx->pos2.first->ch->print();
	printf("=================\n");
	*/

	PG npos1 = opx->pos1;
	PG npos2 = opx->pos2;

	assert(npos1.first != NULL);
	assert(npos1.second != NULL);
	assert(npos2.first != NULL);
	assert(npos2.second != NULL);

	assert(npos1.first->ch == npos1.second->ch);
	assert(npos2.first->ch == npos2.second->ch);

	assert(npos1.first->ch != npos2.first->ch);

	assert(npos1.first->b == npos1.second);
	assert(npos2.first->b == npos2.second);

	npos1.first->b = npos2.second;
	npos2.second->a = npos1.first;
	npos2.first->b = npos1.second;
	npos1.second->a = npos2.first;


	if(npos1.first->ch->type == CIRCULAR && npos2.first->ch->type == CIRCULAR)
	{
		npos2.first->ch->p = NULL;

		remove_empty_chrm(npos2.first->ch);

		gene * q = npos2.second;
		while(true)
		{
			q->ch = npos1.first->ch;
			if(q == npos2.first) break;
			q = q->b;
		}
	}
	else if(npos1.first->ch->type == CIRCULAR)
	{
		assert(npos2.first->ch->type == LINEAR);
		npos1.first->ch->p = NULL;

		remove_empty_chrm(npos1.first->ch);

		gene * q = npos1.second;
		while(true)
		{
			q->ch = npos2.first->ch;
			if(q == npos1.first) break;
			q = q->b;
		}
	}
	else if(npos2.first->ch->type == CIRCULAR)
	{
		assert(npos1.first->ch->type == LINEAR);
		gene * q = npos2.second;
		npos2.first->ch->p = NULL;

		remove_empty_chrm(npos2.first->ch);

		while(true)
		{
			q->ch = npos1.first->ch;
			if(q == npos2.first) break;
			q = q->b;
		}
	}
	else
	{
		assert(npos1.first->ch->type == LINEAR);
		assert(npos2.first->ch->type == LINEAR);

		chrm * ch1 = npos1.first->ch;
		chrm * ch2 = npos2.first->ch;

		gene * q = npos1.first->ch->p;
		while(true)
		{
			q->ch = ch1;
			if(q->b == NULL)
			{
				assert(q->x == 0);
				break;
			}
			q = q->b;
		}

		q = npos2.first->ch->p;
		while(true)
		{
			q->ch = ch2;
			if(q->b == NULL)
			{
				assert(q->x == 0);
				break;
			}
			q = q->b;
		}

	}

	return 0;
}

int genome_base::do_dcj(operation * op)
{
	assert(op->type == DCJ);

	const dcj * opx = dynamic_cast<dcj *>(op);

	PG npos1 = opx->pos1;
	PG npos2 = opx->pos2;

	assert(npos1.first != NULL);
	assert(npos1.second != NULL);
	assert(npos2.first != NULL);
	assert(npos2.second != NULL);

	assert(npos1.first->ch == npos1.second->ch);
	assert(npos2.first->ch == npos2.second->ch);
	assert(npos1.first->b == npos1.second);
	assert(npos2.first->b == npos2.second);

	/*
	if(npos2.first->x == 0 && npos2.second->x == 0)
	{
		// linearization
		assert(npos2.first->a == NULL);
		assert(npos2.second->b == NULL);
		linearization * lop = new linearization(npos1.first, npos1.second, npos2.first, npos2.second);
		do_linearization(lop);
		delete lop;
	}
	*/
	if(npos1.first->ch != npos2.first->ch)
	{
		if(opx->d == false)
		{
			operation * x1 = new reversion(npos2.first);
			do_reversion(x1);
			delete x1;
			npos2 = PG(npos2.second, npos2.first);
		}

		operation * x = new exchange1(npos1, npos2);
		do_exchange1(x);
		delete x;
	}
	else
	{
		bool t = chrm::order(npos1.first, npos2.first);
		if(t == false)
		{
			PG t = npos1;
			npos1 = npos2;
			npos2 = t;
		}

		if(opx->d == true)
		{
			operation * x = new circularization(npos1.second, npos2.first);
			do_circularization(x);
			delete x;
		}
		else
		{
			operation * x = new inversion(npos1.second, npos2.first);
			do_inversion(x);
			delete x;
		}
	}

	return 0;
}

int genome_base::do_circularization(operation * op)
{
	assert(op->type == CIRCULARIZATION);

	const circularization * opx = dynamic_cast<circularization *>(op);

	assert(opx->head != NULL);
	assert(opx->tail != NULL);
	assert(opx->head->a != NULL);
	assert(opx->tail->b != NULL);

	chrm * ch0 = opx->head->ch;

	if(opx->tail->b == opx->head) return 0;
	
	if(opx->head->ch->type == CIRCULAR)
	{
		opx->head->ch->p = opx->tail->b;
	}

	opx->head->a->b = opx->tail->b;
	opx->tail->b->a = opx->head->a;

	opx->head->a = opx->tail;
	opx->tail->b = opx->head;

	chrm * ch = new circular_chrm();

	ch->p = opx->head;

	gene * q = opx->head;
	while(true)
	{
		q->ch = ch;
		if(q == opx->tail) break;
		q = q->b;
	}

	chrms.push_back(ch);

	/*
	if(ch0->p->b->x == 0)
	{
		assert(ch0->type == LINEAR);
		assert(ch0->p->b->b == NULL);
		remove_empty_chrm(ch0);
	}
	*/

	return 0;
}

int genome_base::do_linearization(operation * op)
{
	assert(op->type == LINEARIZATION);

	const linearization * opx = dynamic_cast<linearization *>(op);

	assert(opx->head != NULL);
	assert(opx->tail != NULL);
	assert(opx->head->a != NULL);
	assert(opx->tail->b != NULL);

	assert(opx->head->ch->type == CIRCULAR);
	assert(opx->tail->ch->type == CIRCULAR);

	opx->head->ch->p = NULL;

	remove_empty_chrm(opx->head->ch);

	assert(opx->tail->b == opx->head);

	chrm * ch = new linear_chrm();

	gene * g1 = ch->p;
	gene * g2 = ch->p->b;

	ch->p->b = opx->head;
	opx->head->a = g1;
	g2->a = opx->tail;
	opx->tail->b = g2;
	
	gene * q = ch->p;
	while(true)
	{
		if(q == NULL) break;
		q->ch = ch;
		q = q->b;
	}


	chrms.push_back(ch);

	return 0;
}

int genome_base::do_reversion(operation * op)
{
	assert(op->type == REVERSION);

	const reversion * opx = dynamic_cast<reversion*>(op);

	assert(opx->g != NULL);

	opx->g->ch->reverse();

	return 0;
}

int genome_base::do_inversion(operation * op)
{
	assert(op->type == INVERSION);

	const inversion * opx = dynamic_cast<inversion*>(op);
	assert(opx->head != NULL);
	assert(opx->tail != NULL);
	assert(opx->head->a != NULL);
	assert(opx->tail->b != NULL);

	// update copy
	if(opx->head->a == opx->tail)
	{
		assert(opx->tail->b == opx->head);
		return 0;
	}

	// do operation
	opx->head->a->b = opx->tail;
	opx->tail->b->a = opx->head;

	gene * t = opx->head->a;
	opx->head->a = opx->tail->b;
	opx->tail->b = t;

	gene * q = opx->tail;
	while(true)
	{
		q->reverse();
		if(q == opx->head) break;
		q = q->b;
	}

	return 0;
}

int genome_base::do_insertion(operation * op)
{
	assert(op->type == INSERTION);

	const insertion * opx = dynamic_cast<insertion *>(op);
	
	assert(opx->head != NULL);
	assert(opx->tail != NULL);
	assert(opx->head->a == NULL);
	assert(opx->tail->b == NULL);

	// do operation
	if(opx->dest == NULL)
	{
		assert(opx->head->ch != NULL);
		assert(opx->head->ch->type == CIRCULAR);
		opx->head->a = opx->tail;
		opx->tail->b = opx->head;
		opx->head->ch->p = opx->head;
	}
	else
	{
		assert(opx->dest->a != NULL);

		chrm * ch = opx->dest->ch;

		gene * q = opx->head;
		while(true)
		{
			q->ch = ch;
			if(q == opx->tail) break;
			q = q->b;
		}

		opx->head->a = opx->dest->a;
		opx->dest->a->b = opx->head;
		opx->tail->b = opx->dest;
		opx->dest->a = opx->tail;

	}
	
	return 0;
}

int genome_base::do_deletion(operation * op)
{
	assert(op->type == DELETION);

	const deletion * opx = dynamic_cast<deletion *>(op);

	gene * g1 = opx->head->a;
	gene * g2 = opx->tail->b;
	
	assert(opx->head != NULL);
	assert(opx->tail != NULL);
	assert(opx->head->x != 0);
	assert(opx->tail->x != 0);
	assert(opx->head->a != NULL);
	assert(opx->head->b != NULL);
	assert(opx->tail->a != NULL);
	assert(opx->tail->b != NULL);

	// set the pointer in each chromosome
	opx->tail->b->a = opx->head->a;
	opx->head->a->b = opx->tail->b;

	if(opx->head->ch->type == CIRCULAR)
	{
		if(opx->tail->b == opx->head)
		{
			opx->head->ch->p = NULL;
			remove_empty_chrm(opx->head->ch);
		}
		else
		{
			opx->head->ch->p = opx->tail->b;
		}
	}

	if(opx->b == true)
	{
		gene * q = opx->head;
		while(true)
		{
			gene * next = q->b;
			delete q;
			if(q == opx->tail) break;
			q = next;
		}
	}

	assert(g1 == opx->head->a);
	assert(g2 == opx->tail->b);

	return 0;
}


int genome_base::do_loss(operation * op)
{
	assert(op->type == LOSS);

	const loss * opx = dynamic_cast<loss *>(op);

	operation * x = new deletion(opx->pos, opx->pos);

	operate(x);

	delete x;

	return 0;
}

int genome_base::do_iduplication(operation * op)
{
	assert(op->type == IDUPLICATION);
	
	const iduplication * opx = dynamic_cast<iduplication *>(op);
	assert(opx->dest != NULL);
	assert(opx->head != NULL);
	assert(opx->tail != NULL);

	PG pg = chrm::iduplicate(opx->head, opx->tail);

	operation * x = new insertion(opx->dest, pg.first, pg.second);

	operate(x);

	delete x;

	return 0;
}

int genome_base::do_cduplication(operation * op)
{
	assert(op->type == CDUPLICATION);
	
	const cduplication * opx = dynamic_cast<cduplication *>(op);
	assert(opx->head != NULL);
	assert(opx->tail != NULL);

	chrm * chr = new circular_chrm();
	chrms.push_back(chr);

	PG pg = chrm::duplicate(opx->head, opx->tail);
	operation * x = new insertion(chr->p, pg.first, pg.second);
	operate(x);
	delete x; 

	return 0;
}

int genome_base::do_duplication(operation * op)
{
	assert(op->type == DUPLICATION);
	
	const duplication * opx = dynamic_cast<duplication *>(op);
	assert(opx->dest != NULL);
	assert(opx->head != NULL);
	assert(opx->tail != NULL);

	PG pg = chrm::duplicate(opx->head, opx->tail);

	operation * x = new insertion(opx->dest, pg.first, pg.second);

	operate(x);

	delete x;

	return 0;
}

int genome_base::print() const
{
	for(int i = 0; i < chrms.size(); i++)
	{
		if(chrms.at(i)->p == NULL) continue;
		//if(chrms.at(i)->size() == 0) continue;
		chrms.at(i)->print();
		//chrms.at(i)->reverse_print();
	}
	return 0;
}

int genome_base::print_address() const
{
	for(int i = 0; i < chrms.size(); i++)
	{
		chrms.at(i)->print_address();
	}
	return 0;
}

int genome_base::print_details() const
{
	printf("chromosomes:\n");
	for(int i = 0; i < chrms.size(); i++)
	{
		chrms.at(i)->print();
		//chrms.at(i)->reverse_print();
	}
	printf("\n");
	return 0;
}

/*
pair<int, int> genome_base::pair_encode(int g1, int g2) const
{
	//printf("g1 = %4d, g2 = %4d, alphabet = %4d\n", g1, g2, alphabet_size);

	assert((int)fabs(g1) >= 0 && (int)fabs(g1) <= alphabet_size);
	assert((int)fabs(g2) >= 0 && (int)fabs(g2) <= alphabet_size);

	if((int)fabs(g1) > (int)fabs(g2))
	{
		int t = g1;
		g1 = -1 * g2;
		g2 = -1 * t;
	}
	else if(g1==g2 && g1 < 0)
	{
		g1 = -1 * g1;
		g2 = -1 * g2;
	}

	return pair<int, int>(g1 + alphabet_size, g2 + alphabet_size);
}
*/

gene * genome_base::copy_pointer(const genome_base & gb, gene * p)
{
	assert(p != NULL);
	assert(gb.chrms.size() == chrms.size());

	chrm * ch = p->ch;

	int index = -1;
	for(int i = 0; i < gb.chrms.size(); i++)
	{
		if(gb.chrms.at(i) == ch)
		{
			index = i;
			break;
		}
	}

	assert(index != -1);

	int pos = gb.chrms.at(index)->index(p);

	gene * q = chrms.at(index)->at(pos);

	assert(q->x == p->x);

	return q;
}

int genome_base::remove_empty_chrm(chrm * ch)
{
	assert(ch->p == NULL);
	assert(ch->type == CIRCULAR);

	vector<chrm*>::iterator it;

	for(it = chrms.begin(); it != chrms.end(); it++)
	{
		if(ch == (*it)) break;
	}

	assert(it != chrms.end());

	chrms.erase(it);

	return 0;
}

chrm * genome_base::find_valid_chrm()
{
	for(int i = 0; i < chrms.size(); i++)
	{
		if(chrms.at(i)->p == NULL) continue;
		if(chrms.at(i)->type == CIRCULAR) return chrms.at(i);
		if(chrms.at(i)->p->b->x == 0)
		{
			assert(chrms.at(i)->type == LINEAR);
			assert(chrms.at(i)->p->b->b == NULL);
			continue;
		}

		return chrms.at(i);
	}

	return NULL;
}

dcj * genome_base::copy_dcj(const genome_base & gb, dcj * op)
{
	dcj * opx = new dcj();
	opx->d = op->d;

	opx->pos1.first = copy_pointer(gb, op->pos1.first);
	opx->pos1.second = copy_pointer(gb, op->pos1.second);

	opx->pos2.first = copy_pointer(gb, op->pos2.first);
	opx->pos2.second = copy_pointer(gb, op->pos2.second);

	return opx;
}

reversion * genome_base::copy_reversion(const genome_base & gb, reversion * op)
{
	gene * g = copy_pointer(gb, op->g);
	reversion * opx = new reversion(g);
	return opx;
}
