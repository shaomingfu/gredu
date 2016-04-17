#include "simulator.h"
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <cmath>

using namespace std;

simulator::simulator(config * _conf, genome * _target1, genome * _target2)
	: conf(_conf), source(_conf)
{
	target1 = _target1;
	target2 = _target2;
}

simulator::~simulator()
{
}

int simulator::simulate()
{
	init_source();
	simulate1();
	simulate2();
	return 0;
}

int simulator::init_source()
{
	vector<int> id1;
	vector<int> id2;
	for(int i = 0; i < conf->alphabet_size; i++) id1.push_back(i + 1);
	//for(int i = conf->alphabet_size / 2; i < conf->alphabet_size; i++) id2.push_back(i + 1);

	// TODO
	int r = rand() % 2;
	r = 0;
	if(r==0) source.add_linear_chrm(id1);
	else source.add_circular_chrm(id1);

	// perform random duplications
	for(int k = 0; k < conf->random_duplication_num; k++)
	{
		operation * op = simulate_dup(&source);
		source.operate(op);
	}

	*target1 = source;
	*target2 = source;

	assert(target1->chrms.size() == target2->chrms.size());
	for(int i = 0; i < target1->chrms.size(); i++)
	{
		chrm * sc = target1->chrms.at(i);
		chrm * tc = target2->chrms.at(i);
		gene * sp = sc->p;
		gene * sq = sc->back();
		gene * tp = tc->p;
		gene * tq = tc->back();

		vector<gene*> sv = build_gene_list(PG(sp, sq));
		vector<gene*> tv = build_gene_list(PG(tp, tq));

		assert(sv.size() == tv.size());

		for(int k = 0; k < sv.size(); k++)
		{
			// TODO, includes 0s
			/*
			if(sv.at(k)->x == 0)
			{
				assert(tv.at(k)->x == 0);
				continue;
			}
			*/

			map<gene*, gene*>::iterator sit = s2t.find(sv.at(k));
			map<gene*, gene*>::iterator tit = t2s.find(tv.at(k));
			assert(sit == s2t.end());
			assert(tit == t2s.end());

			s2t.insert(PG(sv.at(k), tv.at(k)));
			t2s.insert(PG(tv.at(k), sv.at(k)));
		}
	}

	return 0;
}

int simulator::simulate1()
{
	for(int i = 0; i < conf->simulation_length; i++)
	{
		operation * op;
		op = simulate_single(target1);

		assert(op != NULL);

		if(op->type == DCJ)
		{
			target1->operate(op);
		}
		else if(op->type == DUPLICATION)
		{
			duplication * opx = dynamic_cast<duplication*>(op);

			PG p(opx->dest->a, opx->dest);
			target1->operate(op);

			vector<gene*> v = build_gene_list(PG(p.first->b, p.second->a));
			for(int k = 0; k < v.size(); k++)
			{
				assert(v.at(k)->x != 0);
				map<gene*, gene*>::iterator it = s2t.find(v.at(k));
				assert(it == s2t.end());
				//s2t.insert(PG(v.at(k), NULL));		// TODO
			}
		}
		else
		{
			assert(1 == 0);
		}
	}
	return 0;
}

int simulator::simulate2()
{
	for(int i = 0; i < conf->simulation_length; i++)
	{
		operation * op;
		op = simulate_single(target2);

		assert(op != NULL);

		if(op->type == DCJ)
		{
			target2->operate(op);
		}
		else if(op->type == DUPLICATION)
		{
			duplication * opx = dynamic_cast<duplication*>(op);

			PG p(opx->dest->a, opx->dest);
			target2->operate(op);

			vector<gene*> v = build_gene_list(PG(p.first->b, p.second->a));
			for(int k = 0; k < v.size(); k++)
			{
				assert(v.at(k)->x != 0);
				map<gene*, gene*>::iterator it = t2s.find(v.at(k));
				assert(it == t2s.end());
				//t2s.insert(PG(v.at(k), NULL));	TODO
			}
		}
		else
		{
			assert(1 == 0);
		}
	}

	return 0;
}

operation* simulator::simulate_dcj(genome * target)
{
	vector<double> inv;
	inv.assign(target->chrms.size(), 0);
	for(int i = 0; i < target->chrms.size(); i++)
	{
		inv.at(i) = (double)(target->chrms.at(i)->count_inversion());
	}

	// INVERSION(DCJ)
	int chrm_index = random(inv);
	if(chrm_index < 0) return NULL;
	assert(target->chrms.at(chrm_index)->size() > 0);

	int start = rand() % target->chrms.at(chrm_index)->size();
	int len = rand() % (target->chrms.at(chrm_index)->size() - start);
	int end = start + len;

	gene * g1 = target->chrms.at(chrm_index)->at(start);
	gene * g2 = target->chrms.at(chrm_index)->at(end);

	return new dcj(PG(g1->a, g1), PG(g2, g2->b), false);
}

operation* simulator::simulate_dup(genome * target)
{
	vector<double> dup;
	dup.assign(target->chrms.size(), 0);
	for(int i = 0; i < target->chrms.size(); i++)
	{
		dup.at(i) = (double)(target->chrms.at(i)->count_duplication(conf->max_duplication_length));
	}

	// DUPLICATION
	int chrm_index = random(dup);
	if(chrm_index < 0) return NULL;

	assert(target->chrms.at(chrm_index)->size() > 0);

	int start, dest;
	int max_length = conf->max_duplication_length;
	int min_length = conf->min_duplication_length;

	if(max_length > target->chrms.at(chrm_index)->size())
	{
		max_length = target->chrms.at(chrm_index)->size();
	}

	int len = min_length + (rand() % (max_length - min_length + 1));

	// TODO, non-overlapping duplication, mismatch with count
	if(target->chrms.at(chrm_index)->type==LINEAR)
	{
		start = rand() % (target->chrms.at(chrm_index)->size() - len + 1);
		dest = rand() % (target->chrms.at(chrm_index)->size() - len + 2);
		dest = (start + len + dest) % (1 + target->chrms.at(chrm_index)->size());
	}
	else
	{
		start = rand() % target->chrms.at(chrm_index)->size();
		dest = rand() % (target->chrms.at(chrm_index)->size() - len + 1);
		dest = (start + len + dest) % (target->chrms.at(chrm_index)->size());
	}

	int end = (start + len - 1) % target->chrms.at(chrm_index)->size();

	return new duplication(target->chrms.at(chrm_index)->at(dest), target->chrms.at(chrm_index)->at(start), target->chrms.at(chrm_index)->at(end));
}

operation* simulator::simulate_single(genome * target)
{
	vector<double> v;
	v.push_back(conf->prob_inversion);
	v.push_back(conf->prob_duplication);
	v.push_back(1.0 - v.at(0) - v.at(1));

	int type = random(v);

	//printf("type = %4d, 0 = %6.3lf, 1 = %6.3lf, 2 = %6.3lf\n", type, v.at(0), v.at(1), v.at(2));

	operation * op;
	if(type==0) op = simulate_dcj(target);
	else if(type==1) op = simulate_dup(target);
	else assert(1 == 0);

	return op;
}
