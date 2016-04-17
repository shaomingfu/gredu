#include "psolver.h"
#include "mygraph.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <algorithm>

using namespace std;

psolver::psolver(genome * g1, genome * g2, int _max_combination)
	: pbase(g1, g2), max_combination(_max_combination)
{
	env = new GRBEnv();

	/*
	print_shared_adjacencies();
	gm1->printi(gi);
	printf("---\n");
	gm2->printi(gi);
	printf("\n");
	*/
}

psolver::~psolver()
{
	delete env;
}

bool psolver::extend()
{
	bool result = false;
	for(int i = 0; i < vsa.size(); i++)
	{
		bool b = extend(vsa[i]);
		if(b == true) result = true;
	}
	return result;
}

bool psolver::simplify1()
{
	bool result = false;
	for(int i = 0; i < cds.size(); i++)
	{
		bool b = simplify1(i);
		if(b == true) result = true;
	}

	return result;
}

bool psolver::simplify2()
{
	bool result = false;
	for(int i = 1; i < gf1.size(); i++)
	{
		if(ex1[i] != NULL) continue;
		assert(ex2[i] == NULL);
		bool b = simplify2(i);
		if(b == true) result = true;
	}
	return result;
}

bool psolver::simplify3()
{
	bool result = false;
	for(int i = 1; i < gf1.size(); i++)
	{
		if(ex1[i] != NULL) continue;
		assert(ex2[i] == NULL);
		bool b1 = simplify3(gf1[i], su1);
		bool b2 = simplify3(gf2[i], su2);
		if(b1 || b2) result = true;
	}
	return result;
}

bool psolver::extend(const shadj & sa)
{
	if(sa.adjacent() == false) return false;

	int f1 = (int)fabs(sa.x1->x);
	int f2 = (int)fabs(sa.x2->x);

	if(f1 == f2) return false;

	if(ex1[f1] != NULL && ex1[f2] != NULL) return false;
	assert(ex2[f1] == NULL || ex2[f2] == NULL);

	if(ex1[f1] != NULL && ex1[f1] != sa.x1) return false;
	if(ex1[f2] != NULL && ex1[f2] != sa.x2) return false;
	if(ex2[f1] != NULL && ex2[f1] != sa.y1) return false;
	if(ex2[f2] != NULL && ex2[f2] != sa.y2) return false;

	if(ex1[f1] == sa.x1 && ex2[f1] == sa.y1)
	{
		if(f2 == 0) return false;
		assert(ex1[f2] == NULL);
		assert(ex2[f2] == NULL);
		ex1[f2] = sa.x2;
		ex2[f2] = sa.y2;
		return true;
	}

	if(ex1[f2] == sa.x2 && ex2[f2] == sa.y2)
	{
		if(f1 == 0) return false;
		assert(ex1[f1] == NULL);
		assert(ex2[f1] == NULL);
		ex1[f1] = sa.x1;
		ex2[f1] = sa.y1;
		return true;
	}

	return false;
}

bool psolver::simplify1(int cdi)
{
	candidate & cd = cds[cdi];

	for(set<int>::const_iterator it = cd.fs.begin(); it != cd.fs.end(); it++)
	{
		int f = *it;
		if(ex1[f] != NULL) return false;
		assert(ex2[f] == NULL);
	}

	vector<int> v;
	v.resize(gf1.size() + 1);
	for(int i = 0; i < cds.size(); i++)
	{
		if(cdi == i) continue;
		vector<int>::iterator it = set_difference(cd.fs.begin(), cd.fs.end(), cds[i].fs.begin(), cds[i].fs.end(), v.begin());
		if(it - v.begin() >= 1) continue;
		if(cds[i].size() > cd.size()) return false;
		assert(cds[i].size() == cd.size());
		if(cds[i].jp == true) return false;
	}

	bool b = testify(cd.xv, cd.yv);
	return b;
}

bool psolver::simplify2(int i)
{
	vector<gene*> xv(1);
	vector<gene*> yv(1);
	for(int j = 0; j < gf1[i].size(); j++)
	{
		xv[0] = gf1[i][j];
		for(int k = 0; k < gf2[i].size(); k++)
		{
			yv[0] = gf2[i][k];
			if(testify(xv, yv) == true) return true;
		}
	}
	return false;
}

bool psolver::simplify3(const vector<gene*> & v, set<gene*> & su)
{
	if(v.size() == 0) return false;

	set<gene*> s;
	for(int i = 0; i < v.size(); i++)
	{
		gene * g = v[i];
		if(vsp[gi[g]].size() > 0) continue;
		if(su.find(g) != su.end()) continue;
		s.insert(g);
	}

	int f = (int)fabs(v[0]->x);

	//printf("family = %3d, num = %3d, non-spanning = %3d\n", f, (int)v.size(), (int)s.size());

	bool b = false;
	for(int i = 0; i < v.size(); i++)
	{
		gene * g = v[i];
		if(su.find(g) != su.end()) continue;
		if(vct[gi[g]].size() >= 1) continue;
		if(s.size() == 0) continue;
		if(s.size() == 1 && s.find(g) != s.end()) continue;
		b = true;

		su.insert(g);
		if(s.find(g) != s.end()) s.erase(g);
		assert(s.size() >= 1);
	}
	
	return b;
}

bool psolver::testify(const vector<gene*> & xv, const vector<gene*> & yv)
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

	if(mi.size() + sav.size() >= max_combination) return false;

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

	//if(nsa == 4) model->write("t.lp");

	model->optimize();
	int code = model->get(GRB_IntAttr_Status);
	int ans = model->get(GRB_DoubleAttr_ObjVal);
	if(code == GRB_CUTOFF) ans = -1;

	printf("seed length = %3d, %4d affected shadjs, %4d affected genes, %4d variables, opt = %3d, nsa = %3d, %s\n",
			(int)xv.size(), (int)sav.size(), (int)mi.size(), (int)vars.size(), ans, nsa, ans <= nsa ? "TRUE" : "FALSE");

	/*
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

	delete model;
	if(ans > nsa) return false;

	for(int i = 0; i < xv.size(); i++)
	{
		gene * x = xv[i];
		gene * y = yv[i];
		int f = (int)fabs(x->x);
		if(f == 0) continue;

		if(ex1[f] != NULL)
		{
			assert(ex1[f] == x);
			assert(ex2[f] == y);
		}

		ex1[f] = x;
		ex2[f] = y;
	}

	return true;
}

set<int> psolver::set_intersection(const set<int> & x, const set<int> & y)
{
	vector<int> z(x.size() + y.size());
	vector<int>::iterator end = std::set_intersection(x.begin(), x.end(), y.begin(), y.end(), z.begin());
	return set<int>(z.begin(), end);
}

set<int> psolver::set_union(const set<int> & x, const set<int> & y)
{
	vector<int> z(x.size() + y.size());
	vector<int>::iterator end = std::set_union(x.begin(), x.end(), y.begin(), y.end(), z.begin());
	return set<int>(z.begin(), end);
}

set<int> psolver::build_span_intersection(const vector<gene*> & v)
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

set<int> psolver::build_span_union(const vector<gene*> & v)
{
	set<int> s;
	for(int i = 0; i < v.size(); i++)
	{
		int g = gi[v[i]];
		s = set_union(s, vsp[g]);
	}
	return s;
}

set<int> psolver::build_contact_union(const vector<gene*> & v)
{
	set<int> s;
	for(int i = 0; i < v.size(); i++)
	{
		int g = gi[v[i]];
		s = set_union(s, vct[g]);
	}
	return s;
}

bool psolver::innocent(const shadj & sa, const vector<gene*> & xv, const vector<gene*> & yv)
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

int psolver::statistic()
{
	int sum = 0;
	for(int i = 1; i < ex1.size(); i++)
	{
		if(ex1[i] != NULL)
		{
			assert(ex2[i] != NULL);
			sum++;
		}
	}

	printf("%d / %d families are fixed\n", sum, (int)ex1.size() - 1);

	return 0;
}

