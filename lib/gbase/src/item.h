#ifndef __ITEM_H__
#define __ITEM_H__

#include "genome_base.h"
#include "adjacency.h"


class item
{
public:
	int type;
	bool valid;
	string name;

public:
	virtual int operate(genome_base & gm) const = 0;
	virtual int reverse_operate(genome_base & gm) const = 0;

	virtual int print() const = 0;

	virtual int find_adjacency(const adjacency & a) = 0;
	virtual int find_gene(gene * a) = 0;
	virtual int make_name() = 0;

	virtual int check_in_active(gene * g) = 0;
	virtual int check_in_active(const adjacency & a) = 0;
	virtual int check_out_active(gene * g) = 0;
	virtual int check_out_active(const adjacency & a) = 0;

};

class gene_item : public item
{
public:
	gene_item(gene * _g);

public:
	gene * g;

public:
	int find_adjacency(const adjacency & a);
	int find_gene(gene * a);
	int operate(genome_base & gm) const;
	int reverse_operate(genome_base & gm) const;
	int print() const;
	int make_name();

	int check_in_active(gene * g);
	int check_in_active(const adjacency & a);
	int check_out_active(gene * g);
	int check_out_active(const adjacency & a);

};

class adj_item : public item
{
public:
	adj_item(const adjacency & _a);

public:
	adjacency a;

public:
	int find_adjacency(const adjacency & a);
	int find_gene(gene * a);
	int operate(genome_base & gm) const;
	int reverse_operate(genome_base & gm) const;
	int print() const;
	int make_name();

	int check_in_active(gene * g);
	int check_in_active(const adjacency & a);
	int check_out_active(gene * g);
	int check_out_active(const adjacency & a);

};

class dcj_item : public item
{
public:
	dcj_item(dcj * opx);
	dcj_item(const adjacency & _a1, const adjacency & _a2, bool _d);
	dcj_item(const adjacency & _a1, const adjacency & _a2, const adjacency & _o1, const adjacency & _o2);
	dcj_item(const dcj_item & di);

public:
	adjacency a1;
	adjacency a2;
	adjacency o1;
	adjacency o2;

public:
	int operate(genome_base & gm) const;
	int reverse_operate(genome_base & gm) const;
	int find_adjacency(const adjacency & a);
	int find_gene(gene * a);
	int print() const;
	int copy(const dcj_item & di);
	int get_map(pair<int, int> & k1, pair<int, int> & k2);
	int get_reverse_map(pair<int, int> & k1, pair<int, int> & k2);
	int make_name();

	//dcj_item reverse();
	int reverse();

	int compare(const dcj_item & x) const;

	int check_in_active(gene * g);
	int check_in_active(const adjacency & a);
	int check_out_active(gene * g);
	int check_out_active(const adjacency & a);

};

class cdup_item : public item
{
public:
	cdup_item();
	cdup_item(const cdup_item & ci);
	cdup_item(const vector<adjacency> & _s, const vector<adjacency> & _t, const adjacency & _c, 
			const vector<gene*>& _sg, const vector<gene*> & _tg);

public:
	vector<adjacency> s;		// source adjacency
	vector<adjacency> t;		// target adjacency
	adjacency c;				// target adjacency

	vector<gene*> sg;			// source genes
	vector<gene*> tg;			// target genes

public:
	int operate(genome_base & gm) const;
	int reverse_operate(genome_base & gm) const;
	int find_adjacency(const adjacency & a);
	int find_gene(gene * a);
	int print() const;
	int copy(const cdup_item & di);

	bool overlap(const cdup_item & x);
	//int compare(const cdup_item & x) const;
	//int rigid_compare(const cdup_item & x) const;
	int make_name();

	int check_in_active(gene * g);
	int check_in_active(const adjacency & a);
	int check_out_active(gene * g);
	int check_out_active(const adjacency & a);

	extremity get_target_boundary_extremity(const extremity & e);
};

class dup_item : public item
{
public:
	dup_item();
	dup_item(const PG & _y, const PG & _z, bool b);

public:
	adjacency y;
	vector<adjacency> s;		// source adjacency
	vector<adjacency> t;		// target adjacency
	adjacency a1, a2;			// two active input adjacencies
	adjacency z;				// active output adjacency

	vector<gene*> sg;			// source genes
	vector<gene*> tg;			// target genes

	int set_x(const PG & x, bool b);
	int set_z(const PG & z, bool b);

public:
	int operate(genome_base & gm) const;
	int reverse_operate(genome_base & gm) const;
	int find_adjacency(const adjacency & a);
	int find_gene(gene * a);
	int print() const;
	int copy(const dup_item & di);

	//bool find_target_gene(gene * g);
	//bool depend(const dup_item & x);
	bool overlap(const dup_item & x);
	int compare(const dup_item & x) const;
	int rigid_compare(const dup_item & x) const;
	int make_name();

	int check_in_active(gene * g);
	int check_in_active(const adjacency & a);
	int check_out_active(gene * g);
	int check_out_active(const adjacency & a);


	pair<item*, item*> split() const;
};

item * copy_item(item * x);

#endif
