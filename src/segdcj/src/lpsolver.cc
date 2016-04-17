#include "lpsolver.h"
#include "gurobi_c++.h"

#include <iostream>
#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <cmath>

lpsolver::lpsolver(config * _conf, genome * _gm1, genome * _gm2)
	: psolver(_conf, _gm1, _gm2) 
{
	env = new GRBEnv();
	model = new GRBModel(*env);
}

lpsolver::~lpsolver()
{
	if(model != NULL) delete model;
	if(env != NULL) delete env;
}

int lpsolver::solve()
{
	// add variables
	build_upper_bounds();

	add_gene_variables(gvars1, gi1.size());
	add_gene_variables(gvars2, gi2.size());

	add_label_variables(lvars1, gi1.size(), vu1);
	add_label_variables(lvars2, gi2.size(), vu2);

	add_bound_variables(xvars1, gi1.size());
	add_bound_variables(xvars2, gi2.size());

	add_edge_variables();

	add_operation_variables(ovars1, dup1.size());
	add_operation_variables(ovars2, dup2.size());

	printf("total (%6d, %6d) variables for genes\n", (int)gvars1.size(), (int)gvars2.size());
	printf("total (%6d, %6d) variables for labels\n", (int)lvars1.size(), (int)lvars2.size());
	printf("total (%6d, %6d) variables for bounds\n", (int)xvars1.size(), (int)xvars2.size());
	printf("total (%6d, %6d) variables for operations\n", (int)ovars1.size(), (int)ovars2.size());
	printf("total (%6d, %6d) variables for gene paris\n", (int)evars.size(), 0);

	// add constraints
	add_gene_label_constraints(gvars1, lvars1, vu1);
	add_gene_label_constraints(gvars2, lvars2, vu2);

	add_adjacency_constraints(gi1, gvars1, lvars1, vu1);
	add_adjacency_constraints(gi2, gvars2, lvars2, vu2);

	add_edge_label_constraints();

	add_bound_constraints(lvars1, xvars1, vu1);
	add_bound_constraints(lvars2, xvars2, vu2);

	add_family_constraints();

	add_gene_edge_constraints();

	add_gene_operation_constraints(gi1, dup1, ovars1, gvars1);
	add_gene_operation_constraints(gi2, dup2, ovars2, gvars2);

	add_operation_constraints(gi1, dup1, ovars1, gvars1, lvars1, vu1);
	add_operation_constraints(gi2, dup2, ovars2, gvars2, lvars2, vu2);

	add_null_gene_constraints(gf1[0], gi1, lvars1);
	add_null_gene_constraints(gf2[0], gi2, lvars2);

	add_inference_constraints();

	set_objective1();

	model->getEnv().set(GRB_DoubleParam_TimeLimit, conf->ilp_timelimit / 3.0);
	model->getEnv().set(GRB_DoubleParam_Heuristics, conf->heuristics);
	model->getEnv().set(GRB_IntParam_MIPFocus, conf->ilp_focus);

	model->update();
	model->optimize();

	set_objective2();

	model->getEnv().set(GRB_DoubleParam_TimeLimit, conf->ilp_timelimit / 1.5);
	model->getEnv().set(GRB_DoubleParam_Heuristics, conf->heuristics);
	model->getEnv().set(GRB_IntParam_MIPFocus, conf->ilp_focus);

	model->update();
	model->optimize();

	// output results
	statistic_solution();
	write_duplicated_segments();
	//print_solution();

	remove_duplicated_genes();
	collect_pairs();

	return 0;
}

int lpsolver::build_upper_bounds()
{
	vu1.assign(ig1.size() * 2, 0);
	vu2.assign(ig2.size() * 2, 0);
	for(int i = 0; i < ig1.size() * 2; i++) vu1[i] = i + 1;
	for(int i = 0; i < ig2.size() * 2; i++) vu2[i] = i + 1 + 2 * ig1.size();
	return 0;
}

int lpsolver::add_gene_variables(vector<GRBVar> & gvars, int n)
{
	gvars.clear();
	for(int i = 0; i < n; i++)
	{
		GRBVar var = model->addVar(0, 1, 0, GRB_BINARY);
		gvars.push_back(var);
	}
	model->update();
	return 0;
}

int lpsolver::add_label_variables(vector<GRBVar> & lvars, int n, const vector<int> & vu)
{
	lvars.clear();
	for(int i = 0; i < n; i++)
	{
		GRBVar var1 = model->addVar(0, vu[2 * i + 0], 0, GRB_INTEGER);
		GRBVar var2 = model->addVar(0, vu[2 * i + 1], 0, GRB_INTEGER);
		lvars.push_back(var1);
		lvars.push_back(var2);
	}
	model->update();
	return 0;
}

int lpsolver::add_bound_variables(vector<GRBVar> & xvars, int n)
{
	xvars.clear();
	for(int i = 0; i < 2 * n; i++)
	{
		GRBVar var = model->addVar(0, 1, 1, GRB_BINARY);
		xvars.push_back(var);
	}
	model->update();
	return 0;
}

int lpsolver::add_edge_variables()
{
	evars.clear();
	for(int i = 0; i < ei.size(); i++)
	{
		GRBVar var = model->addVar(0, 1, 1, GRB_BINARY);
		evars.push_back(var);
	}
	model->update();
	return 0;
}

int lpsolver::add_operation_variables(vector<GRBVar> & ovars, int n)
{
	ovars.clear();
	for(int i = 0; i < n; i++)
	{
		GRBVar var = model->addVar(0, 1, 1, GRB_BINARY);
		ovars.push_back(var);
	}
	model->update();
	return 0;
}

int lpsolver::add_gene_label_constraints(const vector<GRBVar> & gvars, const vector<GRBVar> & lvars, const vector<int> & vu)
{
	for(int i = 0; i < gvars.size(); i++)
	{
		int x = 2 * i + 0;
		int y = 2 * i + 1;
		// if gene is duplicated (gvar = 0), then the two extremieis have the same label
		model->addConstr(lvars[x], GRB_LESS_EQUAL, lvars[y] + gvars[i] * vu[x]);
		model->addConstr(lvars[y], GRB_LESS_EQUAL, lvars[x] + gvars[i] * vu[y]);
		//model->addConstr(lvars[x], GRB_LESS_EQUAL, gvars[i] * vu[x]);
		//model->addConstr(lvars[y], GRB_LESS_EQUAL, gvars[i] * vu[y]);
	}
	return 0;
}

int lpsolver::add_adjacency_constraints(map<gene*, int> & gi, const vector<GRBVar> & gvars, const vector<GRBVar> & lvars, const vector<int> & vu)
{
	map<gene*, int>::const_iterator it;
	for(it = gi.begin(); it != gi.end(); it++)
	{
		gene * x = it->first;
		gene * y = x->b;
		assert(x != NULL);
		if(y == NULL) continue;

		int xi = gi[x];
		int yi = gi[y];
		int xl = xi * 2 + 1;
		int yl = yi * 2 + 0;

		// the two extremities in one adjacency should always have the same label
		model->addConstr(lvars[xl], GRB_EQUAL, lvars[yl]);
		//model->addConstr(lvars[xl], GRB_LESS_EQUAL, lvars[yl] + (2 - gvars[xi] - gvars[yi]) * vu[xl]);
		//model->addConstr(lvars[yl], GRB_LESS_EQUAL, lvars[xl] + (2 - gvars[xi] - gvars[yi]) * vu[yl]);
	}
	return 0;
}

int lpsolver::add_edge_label_constraints()
{
	map<edge_descriptor, int>::const_iterator it;
	for(it = ei.begin(); it != ei.end(); it++)
	{
		int e = it->second;

		int x = source(it->first, pr);
		int y = target(it->first, pr) - ig1.size();
		
		assert(x >= 0 && x < ig1.size());
		assert(y >= 0 && y < ig2.size());

		gene * gx = ig1[x];
		gene * gy = ig2[y];

		assert(gx->x == gy->x || gx->x + gy->x == 0);

		if(gx->x == 0)
		{
			assert(gx->a == NULL || gx->b == NULL);
			assert(gy->a == NULL || gy->b == NULL);

			int xi = (gx->a == NULL) ? 2 * x + 1 : 2 * x + 0;
			int yi = (gy->a == NULL) ? 2 * y + 1 : 2 * y + 0;

			// the two extremities adjacent to the chosen edge should have the same label
			model->addConstr(lvars1[xi], GRB_LESS_EQUAL, lvars2[yi] + (1 - evars[e]) * vu1[xi]);
			model->addConstr(lvars2[yi], GRB_LESS_EQUAL, lvars1[xi] + (1 - evars[e]) * vu2[yi]);
		}
		else
		{
			int xi1 = 2 * x + 0;
			int xi2 = 2 * x + 1;
			int yi1 = (gx->x == gy->x) ? 2 * y + 0 : 2 * y + 1;
			int yi2 = (gx->x == gy->x) ? 2 * y + 1 : 2 * y + 0;

			// the two extremities adjacent to the chosen edge should have the same label
			model->addConstr(lvars1[xi1], GRB_LESS_EQUAL, lvars2[yi1] + (1 - evars[e]) * vu1[xi1]);
			model->addConstr(lvars2[yi1], GRB_LESS_EQUAL, lvars1[xi1] + (1 - evars[e]) * vu2[yi1]);

			model->addConstr(lvars1[xi2], GRB_LESS_EQUAL, lvars2[yi2] + (1 - evars[e]) * vu1[xi2]);
			model->addConstr(lvars2[yi2], GRB_LESS_EQUAL, lvars1[xi2] + (1 - evars[e]) * vu2[yi2]);
		}
	}
	return 0;
}

int lpsolver::add_bound_constraints(const vector<GRBVar> & lvars, const vector<GRBVar> & xvars, const vector<int> & vu)
{
	assert(lvars.size() == xvars.size());
	assert(lvars.size() == vu.size());
	for(int i = 0; i < lvars.size(); i++)
	{
		// guarantee that if xvar = 1 then lvar must be equal to its upper bound
		model->addConstr(vu[i] * xvars[i], GRB_LESS_EQUAL, lvars[i]);
	}
	return 0;
}

int lpsolver::add_family_constraints()
{
	assert(gf1.size() == gf2.size());
	for(int f = 1; f < gf1.size(); f++)
	{
		GRBLinExpr expr1;
		GRBLinExpr expr2;

		for(int i = 0; i < gf1[f].size(); i++)
		{
			gene * g = gf1[f][i];
			int x = gi1[g];
			expr1 += gvars1[x];
		}

		for(int i = 0; i < gf2[f].size(); i++)
		{
			gene * g = gf2[f][i];
			int x = gi2[g];
			expr2 += gvars2[x];
		}

		// for each family, the number of genes that are not duplicated in each genome should be equal
		model->addConstr(expr1, GRB_EQUAL, expr2);

		if(gf1[f].size() >= 1 && gf2[f].size() >= 1)
		{
			model->addConstr(expr1, GRB_GREATER_EQUAL, 1);
			model->addConstr(expr2, GRB_GREATER_EQUAL, 1);
		}
	}
	return 0;
}

int lpsolver::add_gene_edge_constraints()
{
	for(int i = 0; i < gi1.size(); i++)
	{
		out_edge_iterator ei1, ei2;
		GRBLinExpr expr;
		for(tie(ei1, ei2) = out_edges(i, pr); ei1 != ei2; ei1++)
		{
			int e = ei[*ei1];
			expr += evars[e];
		}
		
		// the sum of the edges adjacenct to one gene should be equal to 1 (the gene is not duplicated), 0 (otherwise)
		model->addConstr(expr, GRB_EQUAL, gvars1[i]);
	}

	for(int i = 0; i < gi2.size(); i++)
	{
		out_edge_iterator ei1, ei2;
		GRBLinExpr expr;
		for(tie(ei1, ei2) = out_edges(i + gi1.size(), pr); ei1 != ei2; ei1++)
		{
			int e = ei[*ei1];
			expr += evars[e];
		}

		// the sum of the edges adjacenct to one gene should be equal to 1 (the gene is not duplicated), 0 (otherwise)
		model->addConstr(expr, GRB_EQUAL, gvars2[i]);
	}
	
	return 0;
}

int lpsolver::add_gene_operation_constraints(map<gene*, int> & gi, const vector<PG> & dup, const vector<GRBVar> & ovars, const vector<GRBVar> & gvars)
{
	vector< vector<int> > vv;
	vv.resize(gi.size());
	for(int i = 0; i < dup.size(); i++)
	{
		gene * head = dup[i].first;
		gene * tail = dup[i].second;

		gene * g = head;
		while(true)
		{
			int xi = gi[g];
			vv[xi].push_back(i);
			if(g == tail) break;
			g = g->b;
		}
	}

	for(int i = 0; i < vv.size(); i++)
	{
		GRBLinExpr expr;
		for(int j = 0; j < vv[i].size(); j++)
		{
			expr += ovars[vv[i][j]];
		}

		// each gene can be duplicated at most once
		model->addConstr(expr, GRB_LESS_EQUAL, 1);

		// if a gene is not covered by any duplication, then it must be equal to 1
		model->addConstr(gvars[i], GRB_GREATER_EQUAL, 1 - expr);
	}
	
	return 0;
}

int lpsolver::add_operation_constraints(map<gene*, int> & gi, const vector<PG> & dup, const vector<GRBVar> & ovars,
		const vector<GRBVar> & gvars, const vector<GRBVar> & lvars, const vector<int> & vu)
{
	for(int i = 0; i < dup.size(); i++)
	{
		gene * head = dup[i].first;
		gene * tail = dup[i].second;

		/*
		int ga = gi[head->a];
		int gb = gi[tail->b];
		model->addConstr(gvars[ga], GRB_GREATER_EQUAL, ovars[i]);
		model->addConstr(gvars[gb], GRB_GREATER_EQUAL, ovars[i]);
		*/

		gene * g = head;
		while(true)
		{
			int xi = gi[g];

			// if the operation is chosen, then the target genes should be duplicated
			model->addConstr(gvars[xi], GRB_LESS_EQUAL, 1 - ovars[i]);
			if(g == tail) break;
			g = g->b;
		}

		/*
		int ai = 2 * ga + 1;
		int bi = 2 * gb + 0;
		model->addConstr(lvars[ai], GRB_LESS_EQUAL, lvars[bi] + (1 - ovars[i]) * vu[ai]);
		model->addConstr(lvars[bi], GRB_LESS_EQUAL, lvars[ai] + (1 - ovars[i]) * vu[bi]);
		*/
	}
	return 0;
}

int lpsolver::add_null_gene_constraints(const vector<gene*> & v, map<gene*, int> & gi, const vector<GRBVar> & lvars)
{
	for(int i = 0; i < v.size(); i++)
	{
		gene * g = v[i];
		assert(g->x == 0);
		int x = gi[g];
		int a = x * 2 + 0;
		int b = x * 2 + 1;

		// the two extremities in each null gene should be equal
		model->addConstr(lvars[a], GRB_EQUAL, lvars[b]);
	}
	return 0;
}

int lpsolver::add_inference_constraints()
{
	for(int i = 0; i < vip.size(); i++)
	{
		int x = vip[i].first;
		int y = vip[i].second;
		// the two extremities in each null gene should be equal
		model->addConstr(evars[x], GRB_EQUAL, evars[y]);
	}
	return 0;
}

int lpsolver::set_objective1()
{
	GRBLinExpr expr;
	for(int i = 0; i < gvars1.size(); i++) expr += 0.5 * gvars1.at(i);
	for(int i = 0; i < gvars2.size(); i++) expr += 0.5 * gvars2.at(i);
	for(int i = 0; i < xvars1.size(); i++) expr -= xvars1.at(i);
	for(int i = 0; i < xvars2.size(); i++) expr -= xvars2.at(i);
	for(int i = 0; i < ovars1.size(); i++) expr += conf->sd_weight * ovars1.at(i);
	for(int i = 0; i < ovars2.size(); i++) expr += conf->sd_weight * ovars2.at(i);
	model->setObjective(expr - 1, GRB_MINIMIZE);
	return 0;
}

int lpsolver::set_objective2()
{
	GRBLinExpr expr;
	for(int i = 0; i < gvars1.size(); i++) expr += gvars1.at(i);
	for(int i = 0; i < xvars1.size(); i++) expr -= xvars1.at(i);
	for(int i = 0; i < ovars1.size(); i++) expr += conf->sd_weight * ovars1.at(i);
	for(int i = 0; i < ovars2.size(); i++) expr += conf->sd_weight * ovars2.at(i);
	model->setObjective(expr - 1, GRB_MINIMIZE);
	return 0;
}

int lpsolver::remove_duplicated_genes()
{
	for(int i = 0; i < ovars1.size(); i++)
	{
		double x = ovars1.at(i).get(GRB_DoubleAttr_X);
		if(x <= 0) continue;
		assert(x >= 1);

		operation * op = new deletion(dup1[i].first, dup1[i].second, true);
		gm1->do_deletion(op);
		delete op;
	}

	for(int i = 0; i < ovars2.size(); i++)
	{
		double x = ovars2.at(i).get(GRB_DoubleAttr_X);
		if(x <= 0) continue;
		assert(x >= 1);

		operation * op = new deletion(dup2[i].first, dup2[i].second, true);
		gm2->do_deletion(op);
		delete op;
	}

	return 0;
}

int lpsolver::print()
{
	// print all variables
	printf("gm1:\n");
	gm1->print();
	gm1->statistic();

	printf("operations1: total %6d operations\n", (int)dup1.size());
	/*
	for(int i = 0; i < dup1.size(); i++)
	{
		dup1.at(i)->print();
	}
	*/
	printf("\n");

	printf("gm2:\n");
	gm2->print();
	gm2->statistic();
	printf("operations2: total %6d operations\n", (int)dup2.size());
	/*
	for(int i = 0; i < dup2.size(); i++)
	{
		dup2.at(i)->print();
	}
	*/
	printf("\n");

	return 0;
}

int lpsolver::print_solution()
{
	printf("\ngenome1:\n");
	for(int i = 0; i < gi1.size(); i++)
	{
		gene * g = ig1[i];
		int x = 2 * i + 0;
		int y = 2 * i + 1;
		printf("gene1 %6d: x = %6d, gvar = %1d, upper = (%6d,%6d), lvar = (%6d,%6d), xvar = (%1d,%1d)\n", i, g->x, 
				(int)gvars1[i].get(GRB_DoubleAttr_X),
				vu1[x], vu1[y],
				(int)lvars1[x].get(GRB_DoubleAttr_X),
				(int)lvars1[y].get(GRB_DoubleAttr_X),
				(int)xvars1[x].get(GRB_DoubleAttr_X),
				(int)xvars1[y].get(GRB_DoubleAttr_X));
	}

	printf("\ngenome2:\n");
	for(int i = 0; i < gi2.size(); i++)
	{
		gene * g = ig2[i];
		int x = 2 * i + 0;
		int y = 2 * i + 1;
		printf("gene2 %6d: x = %6d, gvar = %1d, upper = (%6d,%6d), lvar = (%6d,%6d), xvar = (%1d,%1d)\n", i, g->x, 
				(int)gvars2[i].get(GRB_DoubleAttr_X),
				vu2[x], vu2[y],
				(int)lvars2[x].get(GRB_DoubleAttr_X),
				(int)lvars2[y].get(GRB_DoubleAttr_X),
				(int)xvars2[x].get(GRB_DoubleAttr_X),
				(int)xvars2[y].get(GRB_DoubleAttr_X));
	}

	printf("\nedges:\n");
	for(int i = 0; i < ie.size(); i++)
	{
		edge_descriptor e = ie[i];
		int s = source(e, pr);
		int t = target(e, pr) - gi1.size();
		printf("edge  %4d: (%3d,%3d), evar = %1d\n", i, s, t, (int)evars[i].get(GRB_DoubleAttr_X));
	}

	printf("\noperations1:\n");
	for(int i = 0; i < dup1.size(); i++)
	{
		gene * head = dup1[i].first;
		gene * tail = dup1[i].second;

		printf("dup1  %4d:  ovar = %1d, [", i, (int)ovars1[i].get(GRB_DoubleAttr_X));

		gene * g = head;
		while(true)
		{
			int xi = gi1[g];
			printf("(%3d,%3d)", xi, g->x);
			if(g == tail) break;
			g = g->b;
		}
		printf("]\n");
	}

	printf("\noperations2:\n");
	for(int i = 0; i < dup2.size(); i++)
	{
		gene * head = dup2[i].first;
		gene * tail = dup2[i].second;

		printf("dup2  %4d:  ovar = %1d, [", i, (int)ovars2[i].get(GRB_DoubleAttr_X));

		gene * g = head;
		while(true)
		{
			int xi = gi2[g];
			printf("(%3d,%3d)", xi, g->x);
			if(g == tail) break;
			g = g->b;
		}
		printf("]\n");
	}

	return 0;
}

int lpsolver::statistic_solution()
{
	int cycles = 0;
	int genes1 = 0;
	int genes2 = 0;
	map<int, int> m1;
	map<int, int> m2;

	for(int i = 0; i < gi1.size(); i++)
	{
		gene * g = ig1[i];
		if(g->x == 0) continue;
		int x = 2 * i + 0;
		int y = 2 * i + 1;
		if(gvars1[i].get(GRB_DoubleAttr_X) >= 1) genes1++;
		if(xvars1[x].get(GRB_DoubleAttr_X) >= 1) cycles++;
		if(xvars1[y].get(GRB_DoubleAttr_X) >= 1) cycles++;
	}

	for(int i = 0; i < gi2.size(); i++)
	{
		gene * g = ig2[i];
		if(g->x == 0) continue;
		int x = 2 * i + 0;
		int y = 2 * i + 1;
		if(gvars2[i].get(GRB_DoubleAttr_X) >= 1) genes2++;
		if(xvars2[x].get(GRB_DoubleAttr_X) >= 1) cycles++;
		if(xvars2[y].get(GRB_DoubleAttr_X) >= 1) cycles++;
	}

	for(int i = 0; i < dup1.size(); i++)
	{
		int d = distance(dup1[i]);
		if(ovars1[i].get(GRB_DoubleAttr_X) <= 0) continue;
		assert(ovars1[i].get(GRB_DoubleAttr_X) >= 1);
		if(m1.find(d) == m1.end()) m1.insert(PI(d, 1));
		else m1[d]++;
	}

	for(int i = 0; i < dup2.size(); i++)
	{
		int d = distance(dup2[i]);
		if(ovars2[i].get(GRB_DoubleAttr_X) <= 0) continue;
		assert(ovars2[i].get(GRB_DoubleAttr_X) >= 1);
		if(m2.find(d) == m2.end()) m2.insert(PI(d, 1));
		else m2[d]++;
	}

	printf("total %7d cycles in the adjacency graph\n", cycles);
	printf("genome1 total %7d unduplicated genes\n", genes1);
	printf("genome2 total %7d unduplicated genes\n", genes2);

	for(map<int, int>::iterator it = m1.begin(); it != m1.end(); it++)
	{
		printf("genome1 total %7d duplicons of length %3d\n", it->second, it->first + 1);
	}

	for(map<int, int>::iterator it = m2.begin(); it != m2.end(); it++)
	{
		printf("genome2 total %7d duplicons of length %3d\n", it->second, it->first + 1);
	}

	return 0;
}

int lpsolver::write_duplicated_segments()
{
	ofstream fout("segments");
	
	for(int i = 0; i < dup1.size(); i++)
	{
		if(ovars1[i].get(GRB_DoubleAttr_X) <= 0) continue;
		assert(ovars1[i].get(GRB_DoubleAttr_X) >= 1);

		gene * x = dup1[i].first;
		fout<<"G1: ";
		while(true)
		{
			fout<<x->s.c_str()<<" ";	
			if(x == dup1[i].second) break;
			x = x->b;
		}
		fout<<endl;
	}

	for(int i = 0; i < dup2.size(); i++)
	{
		if(ovars2[i].get(GRB_DoubleAttr_X) <= 0) continue;
		assert(ovars2[i].get(GRB_DoubleAttr_X) >= 1);

		gene * x = dup2[i].first;
		fout<<"G2: ";
		while(true)
		{
			fout<<x->s.c_str()<<" ";	
			if(x == dup2[i].second) break;
			x = x->b;
		}
		fout<<endl;
	}

	fout.close();
	return 0;
}

int lpsolver::collect_pairs()
{
	x2y.clear();
	y2x.clear();
	for(int i = 0; i < evars.size(); i++)
	{
		if(evars[i].get(GRB_DoubleAttr_X) <= 0) continue;
		assert(evars[i].get(GRB_DoubleAttr_X) >= 1);
		edge_descriptor e = ie[i];
		int s = source(e, pr);
		int t = target(e, pr);
		gene * x = ig1[s];
		gene * y = ig2[t - gi1.size()];
		x2y.insert(PG(x, y));
		y2x.insert(PG(y, x));
	}
	return 0;
}

int lpsolver::collect_fixed_pairs()
{
	x2y.clear();
	y2x.clear();
	for(map<edge_descriptor, int>::iterator it = ei.begin(); it != ei.end(); it++)
	{
		edge_descriptor e = it->first;
		int s = source(e, pr);
		int t = target(e, pr);
		if(out_degree(s, pr) != 1) continue;
		if(out_degree(t, pr) != 1) continue;
		gene * x = ig1[s];
		gene * y = ig2[t - gi1.size()];
		assert(x2y.find(x) == x2y.end());
		assert(y2x.find(y) == y2x.end());
		x2y.insert(PG(x, y));
		y2x.insert(PG(y, x));
	}
	return 0;
}
