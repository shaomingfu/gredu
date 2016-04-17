#ifndef __OPERATION_H__
#define __OPERATION_H__

#include "gene.h"
#include "common.h"
#include <vector>

using namespace std;

class operation
{
public:
	operation(int _type);

public:
	int type;		// 1: inversion, 2: duplication, 3: loss, more ...

public:
	virtual int print() = 0;
};

typedef pair<operation*, operation*> PO;

class reversion : public operation
{
public:
	reversion();
	reversion(gene * _g);

public:
	gene * g;

public:
	int print();
};

class inversion : public operation
{
public:
	inversion();
	inversion(gene * _head, gene * _tail);

public:
	gene * head;
	gene * tail;

public:
	int print();
};

class linearization: public operation
{
public:
	linearization();
	linearization(gene * _head, gene * _tail);

public:
	gene * head;
	gene * tail;

public:
	int print();
};

class circularization: public operation
{
public:
	circularization();
	circularization(gene * _head, gene * _tail);

public:
	gene * head;
	gene * tail;

public:
	int print();
};

class exchange: public operation
{
public:
	exchange();
	exchange(PG _pos1, PG _ops2);

public:
	PG pos1;
	PG pos2;

public:
	int print();
};


class dcj: public operation
{
public:
	dcj();
	dcj(PG _pos1, PG _ops2, bool _d);

public:
	PG pos1;
	PG pos2;
	bool d;

public:
	int print();
};

class cduplication : public operation
{
public:
	cduplication();
	cduplication(gene * _head, gene * _tail);

public:
	gene * head;
	gene * tail;

public:
	int print();
};

class duplication : public operation
{
public:
	duplication();
	duplication(gene * _dest, gene * _head, gene * _tail);

public:
	gene * dest;
	gene * head;
	gene * tail;

public:
	int print();
};


class iduplication : public operation
{
public:
	iduplication();
	iduplication(gene * _dest, gene * _head, gene * _tail);

public:
	gene * dest;
	gene * head;
	gene * tail;

public:
	int print();
};

class loss : public operation
{
public:
	loss();
	loss(gene * _pos);

public:
	gene * pos;

public:
	int print();

};

class insertion: public operation
{
public:
	insertion();
	insertion(gene * _dest, gene * _head, gene * _tail);

public:
	gene * dest;
	gene * head;
	gene * tail;

public:
	int print();

};

class deletion: public operation
{
public:
	deletion();
	deletion(gene * _head, gene * _tail);
	deletion(gene * _head, gene * _tail, bool _b);

public:
	gene * head;
	gene * tail;
	bool b;

public:
	int print();

};

#endif
