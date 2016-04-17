#include "ilp.h"
#include "mygraph.h"
#include "psolver.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <boost/graph/connected_components.hpp>


using namespace std;

ilp::ilp(config * conf, genome * g1, genome * g2, const MPG & xy, const MPG & yx)
	: pbase(conf, g1, g2, xy, yx)
{
	env = new GRBEnv();
	model = new GRBModel(*env);
	objval = 0;
}

ilp::~ilp()
{
	delete model;
	delete env;
}

int ilp::solve()
{
	/*
	print_shared_adjacencies();
	gm1->printi(gi);
	gm2->printi(gi);
	printf("\n");
	*/

	// variables
	add_gene_variables();
	add_shared_adjacency_variables();

	// constraints
	add_number_constraints();
	add_shared_adjacency_edge_constraints();
	add_coexist_constraints();
	add_inference_constraints();

	set_objective();

	model->getEnv().set(GRB_DoubleParam_TimeLimit, conf->ilp_timelimit);
	model->getEnv().set(GRB_IntParam_MIPFocus, conf->ilp_focus);

	model->update();
	model->optimize();
	objval = model->get(GRB_DoubleAttr_ObjVal);

	collect();
	//analyze();
	update();

	return 0;
}

int ilp::add_gene_variables()
{
	gvars.clear();
	for(int i = 0; i < gi.size(); i++)
	{
		GRBVar var = model->addVar(0, 1, 0, GRB_BINARY);
		gvars.push_back(var);
	}
	model->update();

	printf("total %8d gene variables\n", (int)gvars.size());
	return 0;
}

int ilp::add_shared_adjacency_variables()
{
	avars.clear();
	for(int i = 0; i < vsa.size(); i++)
	{
		GRBVar var = model->addVar(0, 1, 1, GRB_BINARY);
		avars.push_back(var);
	}
	model->update();

	printf("total %8d shared adjacency variables\n", (int)avars.size());
	return 0;
}

int ilp::add_number_constraints()
{

	// keep all the null extremities
	for(int i = 0; i < gf1[0].size(); i++)
	{
		int x = gi[gf1[0][i]];
		model->addConstr(gvars[x], GRB_EQUAL, 1);
	}
	for(int i = 0; i < gf2[0].size(); i++)
	{
		int x = gi[gf2[0][i]];
		model->addConstr(gvars[x], GRB_EQUAL, 1);
	}

	int cnt = gf1[0].size() + gf2[0].size();

	for(int i = 1; i < gf1.size(); i++)
	{
		if(gf1[i].size() == 0 && gf2[i].size() == 0) continue;
		int n = gf1[i].size() < gf2[i].size() ? gf1[i].size() : gf2[i].size();
		if(n >= 1) n = 1;
		GRBLinExpr expr1;
		GRBLinExpr expr2;
		for(int k = 0; k < gf1[i].size(); k++) expr1 += gvars[gi[gf1[i][k]]];
		for(int k = 0; k < gf2[i].size(); k++) expr2 += gvars[gi[gf2[i][k]]];
		model->addConstr(expr1, GRB_EQUAL, expr2);
		model->addConstr(expr1, GRB_GREATER_EQUAL, n);
		model->addConstr(expr2, GRB_GREATER_EQUAL, n);
		cnt += 3;
	}

	printf("total %8d number constraints\n", cnt);
	return 0;
}

int ilp::add_inference_constraints()
{
	int cnt = 0;
	for(int i = 0; i < cds.size(); i++)
	{
		if(cds[i].size() <= 2) continue;
		cnt += cds[i].size() - 2;
		add_inference_constraints(cds[i]);
	}
	printf("total %8d inference constraints\n", cnt);
	return 0;
}

int ilp::add_inference_constraints(const candidate & c)
{
	assert(c.size() >= 3);
	vector<int> s = locate_shared_adjacencies(c);
	if(s.size() <= 1) return 0;
	for(int i = 0; i < s.size() - 1; i++)
	{
		int p = s[i];
		int q = s[i + 1];
		model->addConstr(avars[p], GRB_EQUAL, avars[q]);
	}
	return 0;
}

int ilp::add_shared_adjacency_edge_constraints()
{
	int cnt = 0;
	for(int i = 0; i < vsa.size(); i++)
	{
		shadj & sa = vsa[i];
		assert(gi[sa.x1] < gi[sa.x2]);

		cnt += add_bridge_constraint(sa.x1, sa.x2, i);
		if(sa.direction() == true) cnt += add_bridge_constraint(sa.y1, sa.y2, i);
		else cnt += add_bridge_constraint(sa.y2, sa.y1, i);
	}

	printf("total %8d shared adjacency-gene constraints\n", cnt);
	return 0;
}

int ilp::add_bridge_constraint(gene * p, gene * q, int e)
{
	model->addConstr(gvars[gi[p]], GRB_GREATER_EQUAL, avars[e]);
	model->addConstr(gvars[gi[q]], GRB_GREATER_EQUAL, avars[e]);
	int cnt = 2;
	p = p->b;
	while(p != q)
	{
		model->addConstr(gvars[gi[p]], GRB_LESS_EQUAL, 1 - avars[e]);
		cnt++;
		p = p->b;
	}
	return cnt;
}

int ilp::add_coexist_constraints()
{
	set< pair<int, int> > sp;
	for(int i = 0; i < vct.size(); i++)
	{
		set<int> & s = vct[i];
		set<int>::iterator it1, it2;
		for(it1 = s.begin(); it1 != s.end(); it1++)
		{
			int x = *it1;
			shadj & sax = vsa[x];
			for(it2 = s.begin(); it2 != s.end(); it2++)
			{
				int y = *it2;
				shadj & say = vsa[y];

				if(x >= y) continue;

				bool b = sax.coexist(say);
				if(b == true) continue;

				if(sp.find(pair<int, int>(x, y)) != sp.end()) continue;

				sp.insert(pair<int, int>(x, y));

				model->addConstr(avars[x] + avars[y], GRB_LESS_EQUAL, 1);
			}
		}
	}

	printf("total %8d coexist constraints\n", (int)sp.size());

	return sp.size();
}

int ilp::set_objective()
{
	GRBLinExpr expr;
	for(int i = 0; i < avars.size(); i++) expr += avars.at(i);
	model->setObjective(expr, GRB_MAXIMIZE);
	return 0;
}

int ilp::collect()
{
	saopt.clear();
	for(int k = 0; k < avars.size(); k++)
	{
		if(avars.at(k).get(GRB_DoubleAttr_X) <= 0) continue;
		assert(avars.at(k).get(GRB_DoubleAttr_X) >= 1);
		assert(saopt.find(k) == saopt.end());
		saopt.insert(k);
	}
	return 0;
}

int ilp::update()
{
	// collect fixed pairs
	for(set<int>::iterator it = saopt.begin(); it != saopt.end(); it++)
	{
		int k = *it;
		shadj & sa = vsa[k];
		add_fixed_pair(sa.x1, sa.y1);
		add_fixed_pair(sa.x2, sa.y2);
	}

	// remove non-selected genes
	for(int i = 1; i < gf1.size(); i++)
	{
		for(int j = 0; j < gf1[i].size(); j++)
		{
			gene * g = gf1[i][j];
			int x = gi[g];
			if(gvars.at(x).get(GRB_DoubleAttr_X) >= 1) continue;
			assert(gvars.at(x).get(GRB_DoubleAttr_X) <= 0);
			operation * op = new deletion(g, g, true);
			gm1->do_deletion(op);
			delete op;
		}

		for(int j = 0; j < gf2[i].size(); j++)
		{
			gene * g = gf2[i][j];
			int x = gi[g];
			if(gvars.at(x).get(GRB_DoubleAttr_X) >= 1) continue;
			assert(gvars.at(x).get(GRB_DoubleAttr_X) <= 0);
			operation * op = new deletion(g, g, true);
			gm2->do_deletion(op);
			delete op;
		}
	}

	// map non-fixed genes
	gm1->build_gene_map(gf1);
	gm2->build_gene_map(gf2);
	assert(gf1.size() == gf2.size());

	for(int i = 1; i < gf1.size(); i++)
	{
		assert(gf1[i].size() == gf2[i].size());
		vector<gene*> v1;
		for(int j = 0; j < gf1[i].size(); j++)
		{
			gene * g = gf1[i][j];
			if(is_fixed_gene(g)) continue;
			v1.push_back(g);
		}
		vector<gene*> v2;
		for(int j = 0; j < gf2[i].size(); j++)
		{
			gene * g = gf2[i][j];
			if(is_fixed_gene(g)) continue;
			v2.push_back(g);
		}
		assert(v1.size() == v2.size());

		for(int j = 0; j < v1.size(); j++)
		{
			add_fixed_pair(v1[j], v2[j]);
		}
	}

	return 0;
}

int ilp::analyze()
{
	printf("\n");
	gm1->printv(gi);
	gm2->printv(gi);
	printf("\n");

	psolver p(conf, gm1, gm2, x2y, y2x);
	p.check(saopt);

	p.print_shared_adjacencies(saopt);

	set<int> s;
	for(int k = 0; k < cds.size(); k++)
	{
		vector<int> v = locate_shared_adjacencies(cds[k]);
		for(int i = 0; i < v.size(); i++)
		{
			assert(s.find(v[i]) == s.end());
			s.insert(v[i]);
		}
	}

	set<int> ss = set_intersection(s, saopt);

	print_shared_adjacencies(saopt);

	printf("there are %d adjs in opt, %d in shadj, total %d in shadj\n", (int)saopt.size(), (int)ss.size(), (int)s.size());

	return 0;
}
