#include "ilp6.h"
#include "mygraph.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <boost/graph/connected_components.hpp>


using namespace std;

ilp6::ilp6(const ugraph & _gr, const vector<PII> & _gene_list, const vector<int> & _partners)
	: ilp_base(_gr, _gene_list, _partners)
{
	int n = num_vertices(gr);
	assert(n % 4 == 0);
	max_label_index = 2;
}

ilp6::~ilp6()
{
}

int ilp6::solve()
{

	printf("build neighbor indices ...\n");
	build_neighbor_indices();

	printf("build gene pairs ...\n");
	build_gene_pairs();

	printf("build edge map ...\n");
	build_edge_map();

	set_upper_bounds();

	// variables
	printf("add edge variables ...\n");
	add_edge_variables();

	printf("add label variables ...\n");
	add_label_variables();

	printf("add representative variables ...\n");
	add_representative_variables();

	// constraints
	printf("add degree constraints ...\n");
	add_degree_constraints();

	printf("add connectivity constraints ...\n");
	add_connectivity_constraints();

	printf("add label-representative constraints ...\n");
	add_label_representative_constraints();

	set_objective();

	model->update();
	model->optimize();

	/*
	printf("add extra constraints ...\n");
	add_extra_constraints();

	model->getEnv().set(GRB_IntParam_MIPFocus, 3);
	model->update();
	model->optimize();

	model->getEnv().set(GRB_IntParam_MIPFocus, 1);
	model->update();
	model->optimize();
	*/

	printf("collect results ...\n");
	collect_pairs();
	return 0;
}

int ilp6::set_upper_bounds()
{
	int n = neighbors.size() / 2;
	ub.resize(n);

	for(int i = 0; i < n; i++)
	{
		ub[i] = pow(i, max_label_index);
	}
	return 0;
}

int ilp6::add_edge_variables()
{
	evars.clear();
	char name[1024];
	for(int i = 0; i < gene_pairs.size(); i++)
	{
		sprintf(name, "e%d", i);
		GRBVar var = model->addVar(0, 1, 0, GRB_BINARY, name);
		evars.push_back(var);
	}
	model->update();

	printf("total %8d edge variables\n", (int)evars.size());
	return 0;
}

int ilp6::add_label_variables()
{
	lvars.clear();
	char name[1024];
	assert(neighbors.size() % 2 == 0);
	int n = neighbors.size() / 2;

	for(int i = 0; i < n; i++)
	{
		sprintf(name, "l%d", i);
		GRBVar var = model->addVar(0, ub[i], 0, GRB_CONTINUOUS, name);		// here [0, i] is important
		lvars.push_back(var);
	}
	model->update();
	printf("total %8d label variables\n", (int)lvars.size());
	return 0;
}

int ilp6::add_representative_variables()
{
	rvars.clear();
	char name[1024];
	assert(neighbors.size() % 2 == 0);
	int n = neighbors.size() / 2;

	for(int i = 0; i < n; i++)
	{
		sprintf(name, "r%d", i);
		GRBVar var = model->addVar(0, 1, 1, GRB_BINARY, name);
		rvars.push_back(var);
	}
	model->update();
	printf("total %8d representative variables\n", (int)rvars.size());
	return 0;
}

int ilp6::add_degree_constraints()
{
	char name[1024];
	int cnt = 0;
	out_edge_iterator ei1, ei2;
	for(int i = 0; i < num_vertices(gr); i++)
	{
		GRBLinExpr expr;

		tie(ei1, ei2) = out_edges(i, gr);

		if(eim[*ei1] == -1) assert(distance(ei1, ei2) == 1);
		if(eim[*ei1] == -1) continue;

		for(tie(ei1, ei2) = out_edges(i, gr); ei1 != ei2; ei1++)
		{
			assert(eim[*ei1] >= 0);
			expr += evars.at(eim[*ei1]);
		}
		sprintf(name, "qd%d", i);
		model->addConstr(expr, GRB_EQUAL, 1, name);
		cnt ++;
	}

	printf("total %8d degree constraints\n", cnt);
	return 0;
}

int ilp6::add_extra_constraints()
{
	char name1[1024];
	char name2[1024];
	char name3[1024];

	int cnt = 0;
	int n = num_vertices(gr) / 2;

	out_edge_iterator ei1, ei2;
	out_edge_iterator ei3, ei4;
	for(int i = 0; i < n; i++)
	{
		int ii = nb_indices[i];
		if(ii < 0) continue;
		for(tie(ei1, ei2) = out_edges(i, gr); ei1 != ei2; ei1++)
		{
			int t = target(*ei1, gr);
			assert(t >= n);
			int x = partners[t];
			assert(x >= n);

			GRBLinExpr expr;
			int my = -1;
			for(tie(ei3, ei4) = out_edges(x, gr); ei3 != ei4; ei3++)
			{
				int y = target(*ei3, gr);
				assert(y < n);
				if(i == y) continue;
				
				int yy = nb_indices[y];
				if(ii <= yy) continue;

				if(eim[*ei3] == -1) expr += 1;
				else expr += evars[eim[*ei3]];

				if(yy > my) my = yy;
			}

			if(my != -1)
			{
				double di = pow(ii, max_label_index);
				double dy = pow(my, max_label_index);

				sprintf(name1, "qcx%d.%d", ii, eim[*ei1]);
				sprintf(name2, "qcy%d.%d", ii, eim[*ei1]);
				sprintf(name3, "qcz%d.%d", ii, eim[*ei1]);
				if(eim[*ei1] == -1) 
				{
					model->addConstr(rvars[ii] <= 1 - expr, name1);		
					//model->addConstr(lvars[ii] <= lvars[my] + di * (1 - expr), name2);		
					//model->addConstr(lvars[my] <= lvars[ii] + dy * (1 - expr), name3);		
				}
				else 
				{
					model->addConstr(rvars[ii] <= 2 - evars[eim[*ei1]] - expr, name1);		
					//model->addConstr(lvars[ii] <= lvars[my] + di * (2 - evars[eim[*ei1]] - expr), name2);		
					//model->addConstr(lvars[my] <= lvars[ii] + dy * (2 - evars[eim[*ei1]] - expr), name3);		
				}
				cnt = cnt + 3;
			}
		}
	}
	printf("total %8d extra constraints\n", cnt);
	return 0;
}

int ilp6::add_connectivity_constraints()
{
	char name1[1024];
	char name2[1024];

	int cnt = 0;
	int n = num_vertices(gr) / 2;

	out_edge_iterator ei1, ei2;
	out_edge_iterator ei3, ei4;
	for(int i = 0; i < n; i++)
	{
		int ii = nb_indices[i];
		if(ii < 0) continue;
		for(tie(ei1, ei2) = out_edges(i, gr); ei1 != ei2; ei1++)
		{
			int t = target(*ei1, gr);
			assert(t >= n);
			int x = partners[t];
			assert(x >= n);
			for(tie(ei3, ei4) = out_edges(x, gr); ei3 != ei4; ei3++)
			{
				int y = target(*ei3, gr);
				assert(y < n);
				if(i == y) continue;
				
				int yy = nb_indices[y];
				if(ii >= yy) continue;

				double di = pow(ii, max_label_index);
				double dy = pow(yy, max_label_index);

				//assert(eim[*ei1] >= 0 || eim[*ei3] >= 0);
				sprintf(name1, "qcx%d.%d.%d", ii, eim[*ei1], eim[*ei3]);
				sprintf(name2, "qcy%d.%d.%d", ii, eim[*ei1], eim[*ei3]);
				if(eim[*ei1] == -1 && eim[*ei3] == -1) 
				{
					model->addConstr(lvars[ii] <= lvars[yy], name1);		
					model->addConstr(lvars[yy] <= lvars[ii], name2);		
				}
				else if(eim[*ei1] == -1)
				{
					model->addConstr(lvars[ii] <= lvars[yy] + di * (1 - evars[eim[*ei3]]), name1);		
					model->addConstr(lvars[yy] <= lvars[ii] + dy * (1 - evars[eim[*ei3]]), name2);		
				}
				else if(eim[*ei3] == -1)
				{
					model->addConstr(lvars[ii] <= lvars[yy] + di * (1 - evars[eim[*ei1]]), name1);		
					model->addConstr(lvars[yy] <= lvars[ii] + dy * (1 - evars[eim[*ei1]]), name2);		
				}
				else
				{
					model->addConstr(lvars[ii] <= lvars[yy] + di * (2 - evars[eim[*ei3]] - evars[eim[*ei1]]), name1);		
					model->addConstr(lvars[yy] <= lvars[ii] + dy * (2 - evars[eim[*ei3]] - evars[eim[*ei1]]), name2);		
				}
				cnt = cnt + 2;
			}
		}
	}

	printf("total %8d connectivity constraints\n", cnt);
	return 0;
}

int ilp6::add_label_representative_constraints()
{
	int n = neighbors.size() / 2;
	char name[1024];
	for(int i = 0; i < n; i++)
	{
		sprintf(name, "qlr%d", i);
		model->addConstr(rvars[i] >= lvars[i] - ub[i] + 1, name);		
		model->addConstr(rvars[i] * ub[i] <= lvars[i], name);		
	}
	printf("total %8d label-representative constraints\n", n * 2);
}

int ilp6::set_objective()
{
	GRBLinExpr expr;
	for(int i = 0; i < rvars.size(); i++) expr += rvars.at(i);
	model->setObjective(expr, GRB_MAXIMIZE);
	return 0;
}

int ilp6::collect_pairs()
{
	for(int k = 0; k < gene_pairs.size(); k++)
	{
		if(evars.at(k).get(GRB_DoubleAttr_X) <= 0) continue;
		PII g1 = gene_list[gene_pairs[k].first];
		PII g2 = gene_list[gene_pairs[k].second];
		pairs.insert(PI(g1.first.first, g2.first.first));
		pairs.insert(PI(g1.first.second, g2.first.second));
	}
	return 0;
}

