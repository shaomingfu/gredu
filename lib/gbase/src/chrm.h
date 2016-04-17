#ifndef __CHRM_H__
#define __CHRM_H__

#include "gene.h"
#include <vector>
#include <cstddef>

using namespace std;

class operation;

class chrm
{
public:
	chrm(int _type);
	chrm(const chrm & c);
	virtual ~chrm();

public:
	int type;			// 1 = linear, 2 = circular
	gene * p;

public:
	virtual int size() = 0;
	virtual int print() = 0;
	virtual int blank_print() = 0;
	virtual int print_address() = 0;
	virtual int reverse_print() = 0;

	virtual int regular_print() = 0;
	virtual int regular_print(gene * g) = 0;
	virtual int regular_print(const vector<gene*> & list) = 0;
	virtual int regular_print(const PG & pg) = 0;

	virtual int regular_print(const PG & pg, const PI & pi) = 0;

	virtual int init(const vector<int> & s) = 0;
	virtual int random_init(int len, int alphabet_size) = 0;

	virtual int index(gene * g) = 0;
	virtual gene * at(int index) = 0;
	virtual gene * back() const = 0;
	virtual int destroy() = 0;
	
	virtual int insert_front(int x) = 0;
	virtual int reverse() = 0;

	virtual int count_duplication(int max_length) = 0;
	virtual int count_inversion();
	virtual int count_loss();

	virtual int enumerate_reverse_duplication(chrm * ch, vector<operation*> & v, int max_length) = 0;
	virtual int enumerate_reverse_iduplication(chrm * ch, vector<operation*> & v, int max_length) = 0;
	virtual int enumerate_reverse_duplication(vector<operation*> & v, int max_length) = 0;
	virtual int enumerate_reverse_iduplication(vector<operation*> & v, int max_length) = 0;

	virtual bool is_equal(const chrm * ch) = 0;

	static PG sort(const PG & pg);
	static bool order(gene * x, gene * y);
	static gene * advance(gene * g, int steps);
	static PG make_gene_list(const vector<int> & s, chrm * ch);
	static PG make_gene_list(const vector<int> & s, const vector<string> & ss, chrm * ch);
	static PG duplicate(gene * head, gene * tail);
	static PG iduplicate(gene * head, gene * tail);
	static int print(gene * head, gene * tail);
	static int print(gene * head, int len);
	static int assert_equal(const PG & pg1, const PG& pg2);
	static int reverse(const PG & pg);
	static bool is_equal(const PG & pg1, const PG & pg2);
	static bool is_equal(const PG & pg, const vector<int> & v);
	static bool is_reverse_equal(const PG & pg1, const PG & pg2);
	static bool is_reverse_equal(const PG & pg, const vector<int> & v);
	static bool is_overlap(const PG & pg1, const PG & pg2);
};

#endif
