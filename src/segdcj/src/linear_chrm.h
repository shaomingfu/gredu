#ifndef __LINEAR_CHRM_H__
#define __LINEAR_CHRM_H__

#include "chrm.h"
#include <cstddef>

class linear_chrm : public chrm
{
public:
	linear_chrm();
	linear_chrm(const chrm & lc);
	~linear_chrm();

public:
	int size();
	gene * at(int index);
	int index(gene * g);

	int print();
	int blank_print();
	int print_address();
	int reverse_print();

	int regular_print();
	int regular_print(const vector<gene*> & list);
	int regular_print(const PG & pg);
	int regular_print(const PG & pg, const PI & pi);
	int regular_print(gene * g);

	int insert_front(int x);
	int random_init(int len, int alphabet_size);
	int init(const vector<int> & s);

	gene * back() const;
	int destroy();

	int reverse();

	int count_duplication(int max_length);

	int enumerate_reverse_duplication(chrm * ch, vector<operation*> & v, int max_length);
	int enumerate_reverse_iduplication(chrm * ch, vector<operation*> & v, int max_length);

	int enumerate_reverse_duplication(vector<operation*> & v, int max_length);
	int enumerate_reverse_iduplication(vector<operation*> & v, int max_length);

	bool is_equal(const chrm * ch);
	bool direction(const char * ch);
};

#endif
