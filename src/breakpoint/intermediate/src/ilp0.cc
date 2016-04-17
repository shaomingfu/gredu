#include "ilp0.h"
#include "mygraph.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <boost/graph/connected_components.hpp>


using namespace std;

ilp0::ilp0(config * conf, genome * g1, genome * g2, const MPG & xy, const MPG & yx)
	: pbase0(conf, g1, g2, xy, yx)
{
	env = new GRBEnv();
	model = new GRBModel(*env);
	objval = 0;
}

ilp0::~ilp0()
{
	delete model;
	delete env;
}

int ilp0::solve()
{
	/*
	print_shared_adjacencies();
	gm1->printi(gi);
	gm2->printi(gi);
	printf("\n");
	*/

	// variables
	add_gene_variables();
	add_pair_variables();
	add_shared_adjacency_variables();

	// constraints
	add_number_constraints();
	add_pair_gene_constraints();
	add_fixed_pair_constraints();
	add_pair_coexist_constraints();
	add_shared_adjacency_pair_constraints();
	//add_inference_constraints();

	set_objective();

	model->getEnv().set(GRB_DoubleParam_TimeLimit, conf->ilp_timelimit);
	model->getEnv().set(GRB_IntParam_MIPFocus, conf->ilp_focus);

	model->update();
	model->optimize();
	objval = model->get(GRB_DoubleAttr_ObjVal);

	collect();

	return 0;
}

int ilp0::add_gene_variables()
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

int ilp0::add_pair_variables()
{
	pvars.clear();
	for(int i = 0; i < p2i.size(); i++)
	{
		GRBVar var = model->addVar(0, 1, 0, GRB_BINARY);
		pvars.push_back(var);
	}
	model->update();

	printf("total %8d pair variables\n", (int)pvars.size());
	return 0;
}

int ilp0::add_shared_adjacency_variables()
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

int ilp0::add_number_constraints()
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

int ilp0::add_pair_gene_constraints()
{
	int cnt = 0;
	map<PG, int>::iterator it;
	for(it = p2i.begin(); it != p2i.end(); it++)
	{
		PG p = it->first;
		int k = it->second;

		int gx = gi[p.first];
		int gy = gi[p.second];

		model->addConstr(gvars[gx], GRB_GREATER_EQUAL, pvars[k]);
		model->addConstr(gvars[gy], GRB_GREATER_EQUAL, pvars[k]);
		
		cnt += 2;
	}

	printf("total %8d pair gene constraints\n", cnt);
	return 0;
}

int ilp0::add_fixed_pair_constraints()
{
	int cnt = 0;
	for(MPG::iterator it = x2y.begin(); it != x2y.end(); it++)
	{
		PG pg = PG(it->first, it->second);
		assert(p2i.find(pg) != p2i.end());
		int p = p2i[pg];
		model->addConstr(pvars[p], GRB_EQUAL, 1);
		cnt++;
	}

	printf("total %8d fixed pair constraints\n", cnt);
	return 0;
}

int ilp0::add_pair_coexist_constraints()
{
	int cnt = 0;
	for(int i = 0; i < pct.size(); i++)
	{
		set<int> & s = pct[i];
		if(s.size() <= 1) continue;

		GRBLinExpr expr;
		set<int>::iterator it;
		for(it = s.begin(); it != s.end(); it++) expr += pvars[*it];
		model->addConstr(expr, GRB_LESS_EQUAL, 1);
		cnt++;
	}

	printf("total %8d pair coexist constraints\n", cnt);
	return cnt;
}

int ilp0::add_shared_adjacency_pair_constraints()
{
	int cnt = 0;
	for(int i = 0; i < vsa.size(); i++)
	{
		shadj & sa = vsa[i];
		assert(gi[sa.x1] < gi[sa.x2]);

		int p1 = p2i[PG(sa.x1, sa.y1)];
		int p2 = p2i[PG(sa.x2, sa.y2)];

		model->addConstr(pvars[p1], GRB_GREATER_EQUAL, avars[i]);
		model->addConstr(pvars[p2], GRB_GREATER_EQUAL, avars[i]);

		cnt += add_bridge_constraint(sa.x1, sa.x2, i);
		if(sa.direction() == true) cnt += add_bridge_constraint(sa.y1, sa.y2, i);
		else cnt += add_bridge_constraint(sa.y2, sa.y1, i);
	}

	printf("total %8d shared adjacency pair constraints\n", cnt);
	return 0;
}

int ilp0::add_inference_constraints()
{
	int cnt = 0;
	MPG m = x2y;
	for(int i = 0; i < vsa.size(); i++)
	{
		shadj & sa = vsa[i];
		assert(gi[sa.x1] < gi[sa.x2]);

		int p1 = p2i[PG(sa.x1, sa.y1)];
		int p2 = p2i[PG(sa.x2, sa.y2)];

		vector<gene*> v1;
		vector<gene*> v2;
		v1 = build_gene_list(PG(sa.x1, sa.x2));
		if(sa.direction()) v2 = build_gene_list(PG(sa.y1, sa.y2));
		else v2 = build_gene_list(PG(sa.y2, sa.y1));

		if(v1.size() > 2 || v2.size() > 2) continue;

		GRBLinExpr e = 0;
		for(int k = 1; k < v1.size() - 1; k++) e += gvars[gi[v1[k]]];
		for(int k = 1; k < v2.size() - 1; k++) e += gvars[gi[v2[k]]];

		if(m.find(sa.x1) == m.end() && m.find(sa.y1) == m.end())
		{
			model->addConstr(pvars[p2] - e, GRB_LESS_EQUAL, avars[i]);
			m.insert(PG(sa.x1, sa.y1));
			cnt++;
		}

		if(m.find(sa.x1) != m.end() && m[sa.x1] == sa.y1)
		{
			model->addConstr(pvars[p2] - e, GRB_LESS_EQUAL, avars[i]);
			m.insert(PG(sa.x1, sa.y1));
			cnt++;
		}

		if(m.find(sa.x2) == m.end() && m.find(sa.y2) == m.end())
		{
			model->addConstr(pvars[p1] - e, GRB_LESS_EQUAL, avars[i]);
			m.insert(PG(sa.x2, sa.y2));
			cnt++;
		}

		if(m.find(sa.x2) != m.end() && m[sa.x2] == sa.y2)
		{
			model->addConstr(pvars[p1] - e, GRB_LESS_EQUAL, avars[i]);
			m.insert(PG(sa.x2, sa.y2));
			cnt++;
		}
	}

	printf("total %8d inference constraints\n", cnt);
	return 0;
}

int ilp0::add_bridge_constraint(gene * p, gene * q, int e)
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

int ilp0::set_objective()
{
	GRBLinExpr expr;
	for(int i = 0; i < avars.size(); i++) expr += avars.at(i);
	model->setObjective(expr, GRB_MAXIMIZE);
	return 0;
}

int ilp0::collect()
{
	// collect fixed pairs
	for(int k = 0; k < avars.size(); k++)
	{
		if(avars.at(k).get(GRB_DoubleAttr_X) <= 0) continue;
		assert(avars.at(k).get(GRB_DoubleAttr_X) >= 1);
		shadj & sa = vsa[k];

		assert(gvars.at(gi[sa.x1]).get(GRB_DoubleAttr_X) >= 1);
		assert(gvars.at(gi[sa.x2]).get(GRB_DoubleAttr_X) >= 1);
		assert(gvars.at(gi[sa.y1]).get(GRB_DoubleAttr_X) >= 1);
		assert(gvars.at(gi[sa.y2]).get(GRB_DoubleAttr_X) >= 1);

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
			bool b = is_fixed_gene(g);
			if(b) continue;
			v1.push_back(g);
		}
		vector<gene*> v2;
		for(int j = 0; j < gf2[i].size(); j++)
		{
			gene * g = gf2[i][j];
			bool b = is_fixed_gene(g);
			if(b) continue;
			v2.push_back(g);
		}

		//printf("family %5d, fixed = %2d, total = %2d, v1 = %2d, v2 = %2d\n", i, vf[i], (int)gf1[i].size(), (int)v1.size(), (int)v2.size());

		assert(v1.size() == v2.size());

		for(int j = 0; j < v1.size(); j++)
		{
			add_fixed_pair(v1[j], v2[j]);
		}
	}

	return 0;
}
