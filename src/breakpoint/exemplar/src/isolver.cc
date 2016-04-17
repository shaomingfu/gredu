#include "isolver.h"
#include "mygraph.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <algorithm>

using namespace std;

isolver::isolver(genome * g1, genome * g2)
	: pbase(g1, g2)
{
	env = new GRBEnv();
}

isolver::~isolver()
{
	delete env;
}

bool isolver::testify(const vector<gene*> & xv, const vector<gene*> & yv, vector<int> & sol, const vector<bool> & sab)
{
	GRBModel * model = new GRBModel(*env);

	// construct affected shared-adjacencies 
	set<int> spx = build_span_intersection(xv);
	set<int> spy = build_span_intersection(yv);
	//set<int> spx = build_span_union(xv);
	//set<int> spy = build_span_union(yv);

	set<int> ssa = set_union(spx, spy);

	for(int i = 0; i < xv.size(); i++)
	{
		int f = (int)fabs(xv[i]->x);
		set<int> ctx = build_contact_union(gf1[f]);
		ssa = set_union(ssa, ctx);
	}

	for(int i = 0; i < yv.size(); i++)
	{
		int f = (int)fabs(yv[i]->x);
		set<int> cty = build_contact_union(gf2[f]);
		ssa = set_union(ssa, cty);
	}

	vector<int> sav;
	set<int>::iterator it;
	for(it = ssa.begin(); it != ssa.end(); it++)
	{
		shadj sa = vsa[*it];
		if(innocent(sa, xv, yv) == true) continue;
		if(sab[*it] == false) continue;
		sav.push_back(*it);
	}

	if(sav.size() >= 1 && xv.size() <= 1) return false;

	// construct affected genes
	map<gene*, int> mi;
	for(int i = 0; i < xv.size(); i++)
	{
		int f = (int)fabs(xv[i]->x);
		for(int j = 0; j < gf1[f].size(); j++)
		{
			mi.insert(pair<gene*, int>(gf1[f][j], mi.size()));
		}
	}

	for(int i = 0; i < yv.size(); i++)
	{
		int f = (int)fabs(yv[i]->x);
		for(int j = 0; j < gf2[f].size(); j++)
		{
			mi.insert(pair<gene*, int>(gf2[f][j], mi.size()));
		}
	}

	for(int i = 0; i < sav.size(); i++)
	{
		shadj sa = vsa[sav[i]];
		if(mi.find(sa.x1) == mi.end()) mi.insert(pair<gene*, int>(sa.x1, mi.size()));
		if(mi.find(sa.x2) == mi.end()) mi.insert(pair<gene*, int>(sa.x2, mi.size()));
		if(mi.find(sa.y1) == mi.end()) mi.insert(pair<gene*, int>(sa.y1, mi.size()));
		if(mi.find(sa.y2) == mi.end()) mi.insert(pair<gene*, int>(sa.y2, mi.size()));
	}

	//if(mi.size() + sav.size() >= 500) return false;

	// Variables
	vector<GRBVar> vars;
	for(int i = 0; i < mi.size(); i++)
	{
		GRBVar var = model->addVar(0, 1, 0, GRB_BINARY);
		vars.push_back(var);
	}
	model->update();

	// for each shadj, we add two variables
	for(int i = 0; i < sav.size(); i++)
	{
		GRBVar var1 = model->addVar(0, 1, 1, GRB_BINARY);
		GRBVar var2 = model->addVar(0, 1, 1, GRB_BINARY);
		vars.push_back(var1);
		vars.push_back(var2);
	}
	
	model->update();

	// Constraints
	// each pair of variables for each shadj
	for(int i = 0; i < sav.size(); i++)
	{
		int v1 = mi.size() + i * 2;
		int v2 = mi.size() + i * 2 + 1;
		model->addConstr(vars[v1] + vars[v2], GRB_LESS_EQUAL, 1);
	}

	set<int> sf;
	for(int i = 0; i < xv.size(); i++)
	{
		int f = (int)fabs(xv[i]->x);
		assert(sf.find(f) == sf.end());
		sf.insert(f);
	}

	map<gene*, GRBLinExpr> mg;
	for(int i = 0; i < sav.size(); i++)
	{
		shadj sa = vsa[sav[i]];
		int f1 = (int)fabs(sa.x1->x);
		int f2 = (int)fabs(sa.x2->x);
		if(sf.find(f1) == sf.end() && sf.find(f2) == sf.end()) continue;
		int v1 = mi.size() + i * 2;
		int v2 = mi.size() + i * 2 + 1;
		if(mg.find(sa.x1) == mg.end()) mg.insert(pair<gene*, GRBLinExpr>(sa.x1, vars[v1]));
		else mg[sa.x1] += vars[v1];
		if(mg.find(sa.x2) == mg.end()) mg.insert(pair<gene*, GRBLinExpr>(sa.x2, vars[v2]));
		else mg[sa.x2] += vars[v2];
	}

	map<gene*, GRBLinExpr>::iterator mgi;
	for(mgi = mg.begin(); mgi != mg.end(); mgi++)
	{
		gene * g = mgi->first;
		if(sf.find((int)fabs(g->x)) == sf.end()) model->addConstr(mgi->second, GRB_LESS_EQUAL, 0);
		else model->addConstr(mgi->second, GRB_LESS_EQUAL, 1);
	}

	// do not have the same pair
	for(int i = 0; i < xv.size(); i++)
	{
		if((int)fabs(xv[i]->x) == 0) continue;
		int x = mi[xv[i]];
		int y = mi[yv[i]];
		model->addConstr(vars[x] + vars[y], GRB_LESS_EQUAL, 1);
	}

	// shared adjacencies
	vector<gene*> vx;
	vector<gene*> vy;
	vx.assign(gf1.size(), NULL);
	vy.assign(gf2.size(), NULL);
	for(int i = 0; i < sav.size(); i++)
	{
		shadj sa = vsa[sav[i]];
		if(sa.adjacent() == false) continue;
		int f1 = (int)fabs(sa.x1->x);
		int f2 = (int)fabs(sa.x2->x);
		if(f1 == 0 || f2 == 0) continue;
		if(f1 == f2) continue;

		if( (vx[f2] == NULL || vx[f2] == sa.x2) && (vy[f2] == NULL || vy[f2] == sa.y2) )
		{
			model->addConstr(vars[mi[sa.x1]] + vars[mi[sa.y1]] - 1, GRB_LESS_EQUAL, vars[mi[sa.x2]]);
			model->addConstr(vars[mi[sa.x1]] + vars[mi[sa.y1]] - 1, GRB_LESS_EQUAL, vars[mi[sa.y2]]);
			vx[f2] = sa.x2;
			vy[f2] = sa.y2;
		}

		if( (vx[f1] == NULL || vx[f1] == sa.x1) && (vy[f1] == NULL || vy[f1] == sa.y1) )
		{
			model->addConstr(vars[mi[sa.x2]] + vars[mi[sa.y2]] - 1, GRB_LESS_EQUAL, vars[mi[sa.x1]]);
			model->addConstr(vars[mi[sa.x2]] + vars[mi[sa.y2]] - 1, GRB_LESS_EQUAL, vars[mi[sa.y1]]);
			vx[f1] = sa.x1;
			vy[f1] = sa.y1;
		}
	}

	// shadj --> gene
	for(int i = 0; i < sav.size(); i++)
	{
		shadj sa = vsa[sav[i]];
		model->addConstr(vars[i * 2 + mi.size()] + vars[i * 2 + 1 + mi.size()], GRB_LESS_EQUAL, vars[mi[sa.x1]]);
		model->addConstr(vars[i * 2 + mi.size()] + vars[i * 2 + 1 + mi.size()], GRB_LESS_EQUAL, vars[mi[sa.x2]]);
		model->addConstr(vars[i * 2 + mi.size()] + vars[i * 2 + 1 + mi.size()], GRB_LESS_EQUAL, vars[mi[sa.y1]]);
		model->addConstr(vars[i * 2 + mi.size()] + vars[i * 2 + 1 + mi.size()], GRB_LESS_EQUAL, vars[mi[sa.y2]]);
	}

	map<gene*, int>::iterator mit;
	// conflict using span
	for(mit = mi.begin(); mit != mi.end(); mit++)
	{
		int g = gi[mit->first];
		int x = mit->second;
		for(int i = 0; i < sav.size(); i++)
		{
			if(vsp[g].find(sav[i]) == vsp[g].end()) continue;
			model->addConstr(vars[x] + vars[i * 2 + mi.size()] + vars[i * 2 + 1 + mi.size()], GRB_LESS_EQUAL, 1);
		}
	}

	// gene family size
	map<int, GRBLinExpr> me;
	for(mit = mi.begin(); mit != mi.end(); mit++)
	{
		int f = (int)fabs(mit->first->x);
		if(me.find(f) == me.end())
		{
			GRBLinExpr expr = vars[mit->second];
			me.insert(pair<int, GRBLinExpr>(f, expr));
		}
		else
		{
			me[f] += vars[mit->second];
		}
	}

	map<int, GRBLinExpr>::iterator eit;
	for(eit = me.begin(); eit != me.end(); eit++)
	{
		if(eit->first == 0) continue;
		model->addConstr(eit->second, GRB_LESS_EQUAL, 2);
	}

	// objective function
	GRBLinExpr obj;
	for(int i = 0; i < sav.size(); i++) obj += vars[i * 2 + mi.size()] + vars[i * 2 + 1 + mi.size()];

	int nsa = xv.size() - 1;

	// optimize
	model->getEnv().set(GRB_IntParam_OutputFlag, 0);
	model->setObjective(obj, GRB_MAXIMIZE);
	model->getEnv().set(GRB_IntParam_SolutionLimit, 1);
	model->getEnv().set(GRB_DoubleParam_Cutoff, nsa + 0.5);
	model->update();

	//model->write("t.lp");

	model->optimize();
	int code = model->get(GRB_IntAttr_Status);
	int ans = model->get(GRB_DoubleAttr_ObjVal);
	if(code == GRB_CUTOFF) ans = -1;

	/*
	printf("seed length = %3d, %4d affected shadjs, %4d affected genes, %4d variables, opt = %3d, nsa = %3d, %s\n",
			(int)xv.size(), (int)sav.size(), (int)mi.size(), (int)vars.size(), ans, nsa, ans <= nsa ? "TRUE" : "FALSE");

	printf("segment1: (");
	for(int i = 0; i < xv.size(); i++) printf("%4d[%4d], ", gi[xv[i]], xv[i]->x);
	printf(")\n");

	printf("segment2: (");
	for(int i = 0; i < yv.size(); i++) printf("%4d[%4d], ", gi[yv[i]], yv[i]->x);
	printf(")\n");

	for(int i = 0; i < sav.size(); i++)
	{
		bool f = true;
		if(ans >= 0 && vars[i * 2 + mi.size()].get(GRB_DoubleAttr_X) + vars[i * 2 + 1 + mi.size()].get(GRB_DoubleAttr_X) <= 0) f = false;
		if(f == false) continue;
		string s = vsa[sav[i]].print_string(gi);
		printf("sh-adjs(%4d): %s [%s]\n", sav[i], s.c_str(), f ? "true" : "false");
	}
	printf("\n");
	*/

	sol.clear();
	if(ans > nsa)
	{
		for(int i = 0; i < sav.size(); i++)
		{
			if(vars[i * 2 + mi.size()].get(GRB_DoubleAttr_X) + vars[i * 2 + 1 + mi.size()].get(GRB_DoubleAttr_X) <= 0) continue;
			sol.push_back(sav[i]);

			//string s = vsa[sav[i]].print_string(gi);
			//printf("sh-adjs(%4d): %s\n", sav[i], s.c_str());
		}

		delete model;
		return false;
	}
	else 
	{
		delete model;
		return true;
	}
}

set<int> isolver::set_intersection(const set<int> & x, const set<int> & y)
{
	vector<int> z(x.size() + y.size());
	vector<int>::iterator end = std::set_intersection(x.begin(), x.end(), y.begin(), y.end(), z.begin());
	return set<int>(z.begin(), end);
}

set<int> isolver::set_union(const set<int> & x, const set<int> & y)
{
	vector<int> z(x.size() + y.size());
	vector<int>::iterator end = std::set_union(x.begin(), x.end(), y.begin(), y.end(), z.begin());
	return set<int>(z.begin(), end);
}

set<int> isolver::build_span_intersection(const vector<gene*> & v)
{
	set<int> s;
	vector<int> buf(gi.size());
	if(v.size() == 0) return s;
	int g = gi[v[0]];
	s = vsp[g];
	for(int i = 1; i < v.size(); i++)
	{
		g = gi[v[i]];
		s = set_intersection(s, vsp[g]);
	}
	return s;
}

set<int> isolver::build_span_union(const vector<gene*> & v)
{
	set<int> s;
	for(int i = 0; i < v.size(); i++)
	{
		int g = gi[v[i]];
		s = set_union(s, vsp[g]);
	}
	return s;
}

set<int> isolver::build_contact_union(const vector<gene*> & v)
{
	set<int> s;
	for(int i = 0; i < v.size(); i++)
	{
		int g = gi[v[i]];
		s = set_union(s, vct[g]);
	}
	return s;
}

bool isolver::innocent(const shadj & sa, const vector<gene*> & xv, const vector<gene*> & yv)
{
	assert(xv.size() == yv.size());
	set<int> s;
	for(int i = 0; i < xv.size(); i++)
	{
		int f = (int)fabs(xv[i]->x);
		assert((int)fabs(yv[i]->x) == f);
		assert(s.find(f) == s.end());
		s.insert(f);
	}

	for(int i = 0; i < xv.size(); i++)
	{
		if(xv[i] == sa.x1 && yv[i] == sa.y1 && s.find((int)fabs(sa.x2->x)) == s.end()) return true;
		if(xv[i] == sa.x2 && yv[i] == sa.y2 && s.find((int)fabs(sa.x1->x)) == s.end()) return true;
	}

	return false;
}
