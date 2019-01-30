#include "operation.h"
#include "chrm.h"

#include <cstdio>
#include <cmath>

operation::operation(int _type)
	: type(_type)
{
}

reversion::reversion()
	: operation(REVERSION)
{
	g = NULL;
}

reversion::reversion(gene * _g)
	:operation(REVERSION), g(_g) 
{}

int reversion::print()
{
	printf("reversion: ");
	g->ch->print();
	return 0;
}


inversion::inversion()
	: operation(INVERSION)
{
	head = tail = NULL;
}

inversion::inversion(gene * _head, gene * _tail)
	:operation(INVERSION), head(_head), tail(_tail)
{}

int inversion::print()
{
	printf("inversion: ");
	chrm::print(this->head, this->tail);
	printf("\n");
	return 0;
}

circularization::circularization()
	: operation(CIRCULARIZATION)
{
	head = tail = NULL;
}

circularization::circularization(gene * _head, gene * _tail)
	:operation(CIRCULARIZATION), head(_head), tail(_tail)
{}

int circularization::print()
{
	printf("circularization: (%d, %d)\n", head->x, tail->x);;
	return 0;
}

linearization::linearization()
	: operation(LINEARIZATION)
{
	head = tail = NULL;
}

linearization::linearization(gene * _head, gene * _tail)
	:operation(LINEARIZATION), head(_head), tail(_tail)
{}

int linearization::print()
{
	printf("linearization: (%d, %d)\n", head->x, tail->x);
	return 0;
}

exchange1::exchange1()
	: operation(EXCHANGE)
{
	pos1.first = pos1.second = NULL;
	pos2.first = pos2.second = NULL;
}

exchange1::exchange1(PG _pos1, PG _pos2)
	:operation(EXCHANGE), pos1(_pos1), pos2(_pos2)
{}

int exchange1::print()
{
	printf("exchange1: (%d, %d) + (%d, %d)\n", pos1.first->x, pos1.second->x, pos2.first->x, pos2.second->x);
	return 0;
}


dcj::dcj()
	: operation(DCJ)
{
	pos1.first = pos1.second = NULL;
	pos2.first = pos2.second = NULL;
}

dcj::dcj(PG _pos1, PG _pos2, bool _d)
	:operation(DCJ), pos1(_pos1), pos2(_pos2), d(_d)
{}

int dcj::print()
{
	printf("dcj: {%d} (%d, %d) + (%d, %d)\n", (d==true)?1:0, pos1.first->x, pos1.second->x, pos2.first->x, pos2.second->x);
	return 0;
}

loss::loss()
	: operation(LOSS)
{
	pos = NULL;
}

loss::loss(gene * _pos)
	:operation(LOSS), pos(_pos) 
{}

int loss::print()
{
	printf("loss: ");
	chrm::print(pos, 1);
	printf("\n");
	return 0;
}

insertion::insertion()
	: operation(INSERTION)
{
	dest = NULL;
	head = tail = NULL;
}

insertion::insertion(gene * _dest, gene * _head, gene * _tail)
	:operation(INSERTION), dest(_dest), head(_head), tail(_tail)
{}

int insertion::print()
{
	//printf("insertion       on chromosome %2d, pos   = %3d, size= %3d\n", chrm_index, pos, (int)(genes.size()));
	return 0;
}


deletion::deletion()
	: operation(DELETION)
{
	head = tail = NULL;
	b = true;
}

deletion::deletion(gene * _head, gene * _tail)
	:operation(DELETION), head(_head), tail(_tail)
{
	b = true;
}

deletion::deletion(gene * _head, gene * _tail, bool _b)
	:operation(DELETION), head(_head), tail(_tail)
{
	b = _b;
}

int deletion::print()
{
	//printf("deletion        on chromosome %2d, start = %3d, end = %3d\n", chrm_index, start, end);
	return 0;
}

cduplication::cduplication()
	: operation(CDUPLICATION)
{
	head = tail = NULL;
}

cduplication::cduplication(gene * _head, gene * _tail)
	:operation(CDUPLICATION), head(_head), tail(_tail)
{}

int cduplication::print()
{
	printf(" cduplication: [");
	chrm::print(head, tail);
	printf("]\n");
	return 0;
}

duplication::duplication()
	: operation(DUPLICATION)
{
	dest = NULL;
	head = tail = NULL;
}

duplication::duplication(gene * _dest, gene * _head, gene * _tail)
	:operation(DUPLICATION), dest(_dest), head(_head), tail(_tail)
{}

int duplication::print()
{
	printf(" duplication: [");
	chrm::print(head, tail);
	printf("]\n");
	return 0;

	/*
	printf(" duplication: [%d ", head->a->x);
	chrm::print(head, tail);
	printf(" %d]", tail->b->x);

	printf(" -> [%d ", dest->a->x);
	chrm::print(dest, 1);
	printf("]\n");
	*/
	return 0;
}


iduplication::iduplication()
	: operation(IDUPLICATION)
{
	dest = NULL;
	head = tail = NULL;
}

iduplication::iduplication(gene * _dest, gene * _head, gene * _tail)
	:operation(IDUPLICATION), dest(_dest), head(_head), tail(_tail)
{}

int iduplication::print()
{
	printf("iduplication: [");
	chrm::print(head, tail);
	printf("]\n");
	return 0;

	/*
	printf("iduplication: [%d ", head->a->x);
	chrm::print(head, tail);
	printf(" %d]", tail->b->x);

	printf(" -> [%d ", dest->a->x);
	chrm::print(dest, 1);
	printf("]\n");
	*/

	return 0;
}

