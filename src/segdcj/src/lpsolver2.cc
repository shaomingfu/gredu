#include "lpsolver2.h"
#include "gurobi_c++.h"

#include <iostream>
#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <cmath>

lpsolver2::lpsolver2(config * _conf, const ugraph & _ar, const vector<PI> & _ge,
		const vector<int> & _gf, const vector<int> & _adj,
		const vector< vector<int> > & _vdup, int flag) :
	conf(_conf), ar(_ar), ge(_ge), gf(_gf), adj(_adj), vdup(_vdup)
{
	env = new GRBEnv();
	model = new GRBModel(*env);
	model->getEnv().set(GRB_IntParam_OutputFlag, flag);
	assert(num_vertices(ar) == adj.size());
	assert(ge.size() == gf.size());
}

lpsolver2::~lpsolver2()
{
	if(model != NULL) delete model;
	if(env != NULL) delete env;
}

double lpsolver2::solve()
{
	// data structures
	build_gene_extremities();
	build_edge_map();
	build_upper_bounds();

	// add variables
	add_gene_variables();
	add_label_variables();
	add_bound_variables();
	add_edge_variables();
	add_operation_variables();

	/*
	printf("total %6d variables for genes\n", (int)gvars.size());
	printf("total %6d variables for labels\n", (int)lvars.size());
	printf("total %6d variables for bounds\n", (int)xvars.size());
	printf("total %6d variables for operations\n", (int)ovars.size());
	printf("total %6d variables for gene pairs\n", (int)evars.size());
	*/

	// add constraints
	add_gene_label_constraints();
	add_adjacency_constraints();
	add_edge_label_constraints();
	add_bound_constraints();
	add_family_constraints();
	add_gene_edge_constraints();
	add_gene_operation_constraints();
	add_operation_constraints();

	// set objective function
	set_objective();

	model->update();

	model->getEnv().set(GRB_DoubleParam_TimeLimit, conf->ilp_timelimit);
	model->getEnv().set(GRB_IntParam_MIPFocus, conf->ilp_focus);

	model->optimize();

	// results
	collect_pairs();

	//print_solution();

	return model->get(GRB_DoubleAttr_ObjVal);
}

int lpsolver2::build_gene_extremities()
{
	eg.assign(num_vertices(ar), -1);
	for(int i = 0; i < ge.size(); i++)
	{
		int x = ge[i].first;
		int y = ge[i].second;
		eg[x] = i;
		eg[y] = i;
	}
	return 0;
}

int lpsolver2::build_edge_map()
{
	ei.clear();
	ie.clear();
	edge_iterator ei1, ei2;
	for(tie(ei1, ei2) = edges(ar); ei1 != ei2; ei1++)
	{
		edge_descriptor e = (*ei1);
		int k = ei.size();
		ei.insert(pair<edge_descriptor, int>(e, k));
		ie.insert(pair<int, edge_descriptor>(k, e));
	}
	return 0;
}

int lpsolver2::build_upper_bounds()
{
	vu.assign(num_vertices(ar), 0);
	for(int i = 0; i < num_vertices(ar); i++)
	{
		vu[i] = i + 1;
	}
	return 0;
}

int lpsolver2::add_gene_variables()
{
	gvars.clear();
	for(int i = 0; i < ge.size(); i++)
	{
		GRBVar var = model->addVar(0, 1, 0, GRB_BINARY);
		gvars.push_back(var);
	}
	model->update();
	return 0;
}

int lpsolver2::add_label_variables()
{
	lvars.clear();
	assert(vu.size() == num_vertices(ar));
	for(int i = 0; i < vu.size(); i++)
	{
		GRBVar var = model->addVar(0, vu[i], 0, GRB_INTEGER);
		lvars.push_back(var);
	}
	model->update();
	return 0;
}

int lpsolver2::add_bound_variables()
{
	xvars.clear();
	for(int i = 0; i < vu.size(); i++)
	{
		GRBVar var = model->addVar(0, 1, 1, GRB_BINARY);
		xvars.push_back(var);
	}
	model->update();
	return 0;
}

int lpsolver2::add_edge_variables()
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

int lpsolver2::add_operation_variables()
{
	ovars.clear();
	for(int i = 0; i < vdup.size(); i++)
	{
		GRBVar var = model->addVar(0, 1, conf->sd_weight, GRB_BINARY);
		ovars.push_back(var);
	}
	model->update();
	return 0;
}

int lpsolver2::add_gene_label_constraints()
{
	for(int i = 0; i < gvars.size(); i++)
	{
		int x = ge[i].first;
		int y = ge[i].second;
		// if gene is duplicated (gvar = 0), then the two extremieis have the same label
		model->addConstr(lvars[x], GRB_LESS_EQUAL, lvars[y] + gvars[i] * vu[x]);
		model->addConstr(lvars[y], GRB_LESS_EQUAL, lvars[x] + gvars[i] * vu[y]);
	}
	return 0;
}

int lpsolver2::add_adjacency_constraints()
{
	for(int i = 0; i < adj.size(); i++)
	{
		if(adj[i] == -1) continue;
		int j = adj[i];
		assert(i != j);
		assert(adj[j] == i);
		if(i > j) continue;
		// the two extremities in one adjacency should always have the same label
		model->addConstr(lvars[i], GRB_EQUAL, lvars[j]);
	}
	return 0;
}

int lpsolver2::add_edge_label_constraints()
{
	map<edge_descriptor, int>::iterator it;
	for(it = ei.begin(); it != ei.end(); it++)
	{
		int e = it->second;
		int x = source(it->first, ar);
		int y = target(it->first, ar);
		
		// the two extremities adjacent to the chosen edge should have the same label
		model->addConstr(lvars[x], GRB_LESS_EQUAL, lvars[y] + (1 - evars[e]) * vu[x]);
		model->addConstr(lvars[y], GRB_LESS_EQUAL, lvars[x] + (1 - evars[e]) * vu[y]);
	}
	return 0;
}

int lpsolver2::add_bound_constraints()
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

int lpsolver2::add_family_constraints()
{
	map<int, GRBLinExpr> m;
	for(int i = 0; i < gf.size(); i++)
	{
		if(gf[i] == 0) continue;
		if(m.find(gf[i]) == m.end()) 
		{
			m.insert(pair<int, GRBLinExpr>(gf[i], gvars[i]));
		}
		else
		{
			m[gf[i]] += gvars[i];
		}
	}

	set<int> s;
	map<int, GRBLinExpr>::iterator it;
	for(it = m.begin(); it != m.end(); it++)
	{
		int f = it->first;
		assert(f != 0);
		if(s.find((int)fabs(f)) != s.end()) continue;
		GRBLinExpr expr1 = it->second;
		GRBLinExpr expr2;
		if(m.find(0 - f) != m.end()) expr2 = m[0 - f];

		// for each family, the number of genes that are not duplicated in each genome should be equal
		model->addConstr(expr1, GRB_EQUAL, expr2);

		if(m.find(0 - f) != m.end())
		{
			model->addConstr(expr1, GRB_GREATER_EQUAL, 1);
			model->addConstr(expr2, GRB_GREATER_EQUAL, 1);
		}
	}
	return 0;
}

int lpsolver2::add_gene_edge_constraints()
{
	for(int i = 0; i < num_vertices(ar); i++)
	{
		out_edge_iterator ei1, ei2;
		GRBLinExpr expr;
		for(tie(ei1, ei2) = out_edges(i, ar); ei1 != ei2; ei1++)
		{
			int e = ei[*ei1];
			expr += evars[e];
		}

		if(eg[i] == -1) model->addConstr(expr, GRB_EQUAL, 1);
		else model->addConstr(expr, GRB_EQUAL, gvars[eg[i]]);
	}
	return 0;
}

int lpsolver2::add_gene_operation_constraints()
{
	vector< vector<int> > vv;
	vv.resize(ge.size());
	for(int i = 0; i < vdup.size(); i++)
	{
		for(int j = 0; j < vdup[i].size(); j++)
		{
			int g = vdup[i][j];
			vv[g].push_back(i);
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

int lpsolver2::add_operation_constraints()
{
	for(int i = 0; i < vdup.size(); i++)
	{
		for(int j = 0; j < vdup[i].size(); j++)
		{
			int g = vdup[i][j];
			// if the operation is chosen, then the target genes should be duplicated
			model->addConstr(gvars[g], GRB_LESS_EQUAL, 1 - ovars[i]);
		}
	}
	return 0;
}

int lpsolver2::add_consistent_constraints()
{
	map<edge_descriptor, int>::iterator it;
	for(it = ei.begin(); it != ei.end(); it++)
	{
		int e = it->second;
		int s = source(it->first, ar);
		int t = target(it->first, ar);
		int x = eg[s];
		int y = eg[t];
		if(x == -1) assert(y == -1);
		if(y == -1) assert(x == -1);
		if(x == -1) continue;

		int a = -1;
		int b = -1;
		if(ge[x].first == s)
		{
			assert(ge[y].first == t);
			a = ge[x].second;
			b = ge[y].second;
		}
		else if(ge[x].second == s)
		{
			assert(ge[y].second == t);
			a = ge[x].first;
			b = ge[y].first;
		}
		else assert(false);

		pair<edge_descriptor, int> p = edge(a, b, ar);
		assert(p.second == true);

		int c = ei[p.first];
		assert(e != c);
		if(e > c) continue;

		model->addConstr(evars[e], GRB_EQUAL, evars[c]);
	}
	return 0;
}

int lpsolver2::set_objective()
{
	GRBLinExpr expr;
	for(int i = 0; i < ovars.size(); i++) expr += conf->sd_weight * ovars.at(i);
	for(int i = 0; i < gvars.size(); i++) expr += 0.5 * gvars.at(i);
	for(int i = 0; i < xvars.size(); i++) expr -= xvars.at(i);
	model->setObjective(expr - 1, GRB_MINIMIZE);
	return 0;
}

int lpsolver2::print_solution()
{
	for(int i = 0; i < num_vertices(ar); i++)
	{
		int g = eg[i];
		int dup = (g == -1) ? 0 : 1 - gvars[g].get(GRB_DoubleAttr_X);
		printf("vertex %5d, label = %5d, upper bound = %5d, representative = %1d, gene = %5d, duplicated = %1d\n",
				i, (int)lvars[i].get(GRB_DoubleAttr_X), vu[i], (int)xvars[i].get(GRB_DoubleAttr_X), eg[i], dup);
	}
	return 0;
}

int lpsolver2::get_num_cycles()
{
	int c = 0;
	for(int i = 0; i < xvars.size(); i++)
	{
		double x = xvars[i].get(GRB_DoubleAttr_X);
		if(x <= 0) continue;
		assert(x >= 1);
		c++;
	}
	return c;
}

int lpsolver2::collect_pairs()
{
	x2y.clear();
	y2x.clear();
	for(int i = 0; i < evars.size(); i++)
	{
		if(evars[i].get(GRB_DoubleAttr_X) <= 0) continue;
		assert(evars[i].get(GRB_DoubleAttr_X) >= 1);
		edge_descriptor e = ie[i];
		int s = source(e, ar);
		int t = target(e, ar);
		int x = eg[s];
		int y = eg[t];
		x2y.insert(PI(x, y));
		y2x.insert(PI(y, x));
	}
	return 0;
}
