#ifndef __GENOME_BASE_H__
#define __GENOME_BASE_H__

#include "chrm.h"
#include "operation.h"

#include <vector>
#include <map>

using namespace std;

class genome_base
{
public:
	genome_base(int _alphabet_size);
	genome_base(const genome_base & g);
	genome_base & operator= (const genome_base & g);
	virtual ~genome_base();

public:
	int alphabet_size;
	vector<chrm*> chrms;

public:
	// add chromosome
	int add_linear_chrm(const vector<int> & s);
	int add_circular_chrm(const vector<int> & s);
	int add_linear_chrm(const vector<int> & s, const vector<string> & ss);
	int add_circular_chrm(const vector<int> & s, const vector<string> & ss);

	// do the operation
	int operate(const PO & po);
	int operate(operation * op);
	int do_inversion(operation * op);
	int do_insertion(operation * op);
	int do_deletion(operation * op);
	int do_cduplication(operation * op);
	int do_duplication(operation * op);
	int do_iduplication(operation * op);
	int do_loss(operation * op);
	int do_circularization(operation * op);
	int do_exchange(operation * op);
	int do_dcj(operation * op);
	int do_reversion(operation * op);
	int do_linearization(operation * op);

	// calculate the reverse operation
	PO reverse(operation * op) const;
	PO reverse_dcj(operation * op) const;
	PO reverse_reversion(operation * op) const;
	PO reverse_circularization(operation * op) const;
	PO reverse_linearization(operation * op) const;

	// print out all chrms
	int print() const;
	int print_address() const;
	int print_details() const;

	chrm * find_valid_chrm();

	gene * copy_pointer(const genome_base & gb, gene * p);
	dcj * copy_dcj(const genome_base & gb, dcj * op);
	reversion * copy_reversion(const genome_base & gb, reversion * op);

public:
	int remove_empty_chrm(chrm * ch);
	pair<int, int> pair_encode(int g1, int g2) const;
};

#endif
