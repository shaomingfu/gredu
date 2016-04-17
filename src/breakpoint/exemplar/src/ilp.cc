#include "ilp.h"
#include "mygraph.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <boost/graph/connected_components.hpp>


using namespace std;

ilp::ilp(genome * g1, genome * g2)
	: spliter(g1, g2)
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
	add_inference_constraints();

	set_objective();

	model->update();
	model->optimize();
	objval = model->get(GRB_DoubleAttr_ObjVal);

	collect();
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
	int cnt = 0;
	for(int i = 1; i < gf1.size(); i++)
	{
		if(gf1[i].size() == 0 && gf2[i].size() == 0) continue;
		int n = gf1[i].size() < gf2[i].size() ? gf1[i].size() : gf2[i].size();
		if(n >= 1) n = 1;
		GRBLinExpr expr1;
		GRBLinExpr expr2;
		for(int k = 0; k < gf1[i].size(); k++) expr1 += gvars[gi[gf1[i][k]]];
		for(int k = 0; k < gf2[i].size(); k++) expr2 += gvars[gi[gf2[i][k]]];
		model->addConstr(expr1, GRB_EQUAL, n);
		model->addConstr(expr2, GRB_EQUAL, n);
		cnt += 2;
	}

	printf("total %8d number constraints\n", cnt);
	return 0;
}

int ilp::add_inference_constraints()
{
	int cnt = 0;
	for(int i = 0; i < cps.size(); i++)
	{
		cnt += cps[i].size() - 1;
		add_inference_constraints(cps[i]);
	}
	printf("total %8d inference constraints\n", cnt);
	return 0;
}

int ilp::add_inference_constraints(const candidate & c)
{
	for(int i = 0; i < c.size() - 1; i++)
	{
		int f1 = (int)fabs(c.xv[i]->x);
		int f2 = (int)fabs(c.xv[i + 1]->x);

		if(f1 == 0 || f2 == 0) continue;
		assert(f1 != f2);

		int s1 = gi[c.xv[i]];
		int s2 = gi[c.xv[i + 1]];
		int t1 = gi[c.yv[i]];
		int t2 = gi[c.yv[i + 1]];

		model->addConstr(gvars[s2] + gvars[t2] - 1, GRB_LESS_EQUAL, gvars[s1]);
		model->addConstr(gvars[s2] + gvars[t2] - 1, GRB_LESS_EQUAL, gvars[t1]);

		model->addConstr(gvars[s1] + gvars[t1] - 1, GRB_LESS_EQUAL, gvars[s2]);
		model->addConstr(gvars[s1] + gvars[t1] - 1, GRB_LESS_EQUAL, gvars[t2]);
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

int ilp::set_objective()
{
	GRBLinExpr expr;
	for(int i = 0; i < avars.size(); i++) expr += avars.at(i);
	model->setObjective(expr, GRB_MAXIMIZE);
	return 0;
}

int ilp::collect()
{
	set<gene*> s;
	for(int k = 0; k < gvars.size(); k++)
	{
		if(gvars.at(k).get(GRB_DoubleAttr_X) <= 0) continue;
		assert(gvars.at(k).get(GRB_DoubleAttr_X) >= 1);
		s.insert(ig[k]);
	}

	ex1.resize(gf1.size(), NULL);
	ex2.resize(gf2.size(), NULL);

	for(int f = 1; f < gf1.size(); f++)
	{
		gene * x = NULL;
		for(int k = 0; k < gf1[f].size(); k++)
		{
			gene * g = gf1[f][k];
			if(s.find(g) != s.end())
			{
				assert(x == NULL);
				x = g;
			}
		}

		gene * y = NULL;
		for(int k = 0; k < gf2[f].size(); k++)
		{
			gene * g = gf2[f][k];
			if(s.find(g) != s.end())
			{
				assert(y == NULL);
				y = g;
			}
		}

		if(x == NULL) assert(y == NULL);
		if(y == NULL) assert(x == NULL);
		
		ex1[f] = x;
		ex2[f] = y;
	}

	return 0;
}
