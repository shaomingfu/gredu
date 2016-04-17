#include "psolver.h"
#include "verifier.h"
#include "mygraph.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <algorithm>

using namespace std;

psolver::psolver(config * conf, genome * g1, genome * g2, const MPG & xy, const MPG & yx)
	: pbase(conf, g1, g2, xy, yx)
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
	for(int f = 1; f < gf1.size(); f++)
	{
		bool b = simplify2(f);
		if(b == true) result = true;
	}
	return result;
}

bool psolver::extend(const shadj & sa)
{
	if(sa.adjacent() == false) return false;

	if(is_fixed_pair(sa.x1, sa.y1) == true)
	{
		if(is_fixed_gene(sa.x2)) return false;
		if(is_fixed_gene(sa.y2)) return false;
		add_fixed_pair(sa.x2, sa.y2);
		return true;
	}

	if(is_fixed_pair(sa.x2, sa.y2) == true)
	{
		if(is_fixed_gene(sa.x1)) return false;
		if(is_fixed_gene(sa.y1)) return false;
		add_fixed_pair(sa.x1, sa.y1);
		return true;
	}

	return false;
}

bool psolver::simplify1(int cdi)
{
	candidate & cd = cds[cdi];

	for(int k = 0; k < cd.xv.size(); k++)
	{
		if(is_fixed_gene(cd.xv[k])) return false;
		if(is_fixed_gene(cd.yv[k])) return false;
	}

	verifier vr(cd.xv, cd.yv, cd.d, this, env, conf->max_combination);
	bool b = vr.verify();

	return b;
}

bool psolver::simplify2(int f)
{
	vector<gene*> xv(1);
	vector<gene*> yv(1);
	for(int j = 0; j < gf1[f].size(); j++)
	{
		xv[0] = gf1[f][j];
		if(is_fixed_gene(xv[0]) == true) continue;
		for(int k = 0; k < gf2[f].size(); k++)
		{
			yv[0] = gf2[f][k];
			if(is_fixed_gene(yv[0]) == true) continue;

			verifier vr(xv, yv, true, this, env, conf->max_combination);
			bool b = vr.verify();
			if(b == true) return true;
		}
	}
	return false;
}

bool psolver::testify(const vector<gene*> & xv, const vector<gene*> & yv)
{

	// construct affected shared-adjacencies 
	set<int> spx = build_span_intersection(xv);
	set<int> spy = build_span_intersection(yv);

	set<int> ctx = build_contact_union(xv);
	set<int> cty = build_contact_union(yv);

	set<int> ssa;
	ssa = set_union(ssa, spx);
	ssa = set_union(ssa, spy);
	ssa = set_union(ssa, ctx);
	ssa = set_union(ssa, cty);

	vector<int> v0 = locate_shared_adjacencies(xv, yv);

	// construct affected shared adjacencies
	map<int, int> sam2;		// affected shared adjacencies
	vector<int> sav2;
	set<int>::iterator it;
	for(it = ssa.begin(); it != ssa.end(); it++)
	{
		shadj & sa = vsa[*it];
		if(innocent(sa, xv, yv) == true) continue;
		sam2.insert(PI(*it, sav2.size()));
		sav2.push_back(*it);
	}

	if(sav2.size() >= 1 && xv.size() <= 1) return false;

	// construct affected genes
	map<gene*, int> mi;
	for(int i = 0; i < xv.size(); i++)
	{
		mi.insert(pair<gene*, int>(xv[i], mi.size()));
		mi.insert(pair<gene*, int>(yv[i], mi.size()));
	}

	for(int i = 0; i < sav2.size(); i++)
	{
		shadj & sa = vsa[sav2[i]];
		vector<gene*> v = sa.get_genes();
		for(int k = 0; k < v.size(); k++)
		{
			mi.insert(pair<gene*, int>(v[k], mi.size()));
		}
	}

	// construct innocent shared adjacencies
	map<int, int> sam1;		// affected shared adjacencies
	vector<int> sav1;
	for(it = ssa.begin(); it != ssa.end(); it++)
	{
		shadj & sa = vsa[*it];
		if(innocent(sa, xv, yv) == false) continue;
		vector<gene*> v = sa.get_genes();
		bool b = true;
		for(int k = 0; k < v.size(); k++)
		{
			gene * g = v[k];
			if(mi.find(g) == mi.end()) b = false;
			if(b == false) break;
		}
		if(b == false) continue;

		sam1.insert(PI(*it, sav1.size()));
		sav1.push_back(*it);
	}

	for(int k = 0; k < v0.size(); k++)
	{
		int x = v0[k];
		assert(sam1.find(x) != sam1.end());
		assert(sam2.find(x) == sam2.end());
	}

	if(mi.size() + sav2.size() >= conf->max_combination) return false;
	if(mi.size() + sav1.size() >= conf->max_combination) return false;

	// PART1, compute the optimal innocent shared adjacencies

	// Model
	GRBModel * model1 = new GRBModel(*env);

	// Variables
	// for each involved gene
	vector<GRBVar> vars1;
	for(int i = 0; i < mi.size(); i++)
	{
		GRBVar var = model1->addVar(0, 1, 0, GRB_BINARY);
		vars1.push_back(var);
	}
	model1->update();

	// for each involved innocent shared adjacency
	for(int i = 0; i < sav1.size(); i++)
	{
		GRBVar var = model1->addVar(0, 1, 1, GRB_BINARY);
		vars1.push_back(var);
	}
	
	model1->update();

	// Constraints
	// given shared segments are automatically selected
	for(int k = 0; k < v0.size(); k++)
	{
		int x = v0[k];
		assert(sam1.find(x) != sam1.end());
		int xi = sam1[x];
		model1->addConstr(vars1[xi + mi.size()], GRB_EQUAL, 1);
	}

	// shadj --> boundary gene
	for(int i = 0; i < sav1.size(); i++)
	{
		shadj sa = vsa[sav1[i]];
		model1->addConstr(vars1[i + mi.size()], GRB_LESS_EQUAL, vars1[mi[sa.x1]]);
		model1->addConstr(vars1[i + mi.size()], GRB_LESS_EQUAL, vars1[mi[sa.x2]]);
		model1->addConstr(vars1[i + mi.size()], GRB_LESS_EQUAL, vars1[mi[sa.y1]]);
		model1->addConstr(vars1[i + mi.size()], GRB_LESS_EQUAL, vars1[mi[sa.y2]]);
	}

	map<gene*, int>::iterator mit;

	// shadj --> inbetween genes
	set<int> sas1(sav1.begin(), sav1.end());
	for(mit = mi.begin(); mit != mi.end(); mit++)
	{
		int g = gi[mit->first];
		int x = mit->second;
		set<int> s = set_intersection(vsp[g], sas1);

		for(set<int>::iterator it = s.begin(); it != s.end(); it++)
		{
			assert(sam1.find(*it) != sam1.end());
			int k = sam1[*it];
			model1->addConstr(vars1[x] + vars1[k + mi.size()], GRB_LESS_EQUAL, 1);
		}
	}

	// shadj coexist constraints
	set<PI> sp;
	for(mit = mi.begin(); mit != mi.end(); mit++)
	{
		int g = gi[mit->first];
		set<int> s = set_intersection(vct[g], sas1);
		set<int>::iterator it1, it2;
		for(it1 = s.begin(); it1 != s.end(); it1++)
		{
			int x = *it1;
			shadj & sax = vsa[x];
			assert(sam1.find(x) != sam1.end());
			int xi = sam1[x];
			for(it2 = s.begin(); it2 != s.end(); it2++)
			{
				int y = *it2;
				shadj & say = vsa[y];
				assert(sam1.find(y) != sam1.end());
				int yi = sam1[y];

				if(x >= y) continue;
				bool b = sax.coexist(say);
				if(b == true) continue;

				if(sp.find(pair<int, int>(x, y)) != sp.end()) continue;
				sp.insert(pair<int, int>(x, y));

				model1->addConstr(vars1[mi.size() + xi] + vars1[mi.size() + yi], GRB_LESS_EQUAL, 1);
			}
		}
	}

	// objective function
	GRBLinExpr obj1;
	for(int i = 0; i < sav1.size(); i++) obj1 += vars1[i + mi.size()];

	// optimize
	model1->getEnv().set(GRB_IntParam_OutputFlag, 0);
	model1->setObjective(obj1, GRB_MAXIMIZE);
	model1->update();

	model1->optimize();
	int ans1 = model1->get(GRB_DoubleAttr_ObjVal);

	delete model1;

	printf("given %3d shared adjacencies, optimal value is %3d\n", (int)v0.size(), ans1);

	// PART 2, compute the optimal affected shared adjacencies
	GRBModel * model2 = new GRBModel(*env);


	// Variables
	// for each involved gene
	vector<GRBVar> vars2;
	for(int i = 0; i < mi.size(); i++)
	{
		GRBVar var = model2->addVar(0, 1, 0, GRB_BINARY);
		vars2.push_back(var);
	}
	model2->update();

	// for each involved shared adjacency
	for(int i = 0; i < sav2.size(); i++)
	{
		GRBVar var = model2->addVar(0, 1, 1, GRB_BINARY);
		vars2.push_back(var);
	}
	
	model2->update();

	// Constraints
	// shadj --> boundary gene
	for(int i = 0; i < sav2.size(); i++)
	{
		shadj & sa = vsa[sav2[i]];
		assert(mi.find(sa.x1) != mi.end());
		assert(mi.find(sa.x2) != mi.end());
		assert(mi.find(sa.y1) != mi.end());
		assert(mi.find(sa.y2) != mi.end());

		model2->addConstr(vars2[i + mi.size()], GRB_LESS_EQUAL, vars2[mi[sa.x1]]);
		model2->addConstr(vars2[i + mi.size()], GRB_LESS_EQUAL, vars2[mi[sa.x2]]);
		model2->addConstr(vars2[i + mi.size()], GRB_LESS_EQUAL, vars2[mi[sa.y1]]);
		model2->addConstr(vars2[i + mi.size()], GRB_LESS_EQUAL, vars2[mi[sa.y2]]);
	}


	// shadj --> inbetween genes
	set<int> sas2(sav2.begin(), sav2.end());
	for(mit = mi.begin(); mit != mi.end(); mit++)
	{
		int g = gi[mit->first];
		int x = mit->second;
		set<int> s = set_intersection(vsp[g], sas2);

		for(set<int>::iterator it = s.begin(); it != s.end(); it++)
		{
			assert(sam2.find(*it) != sam2.end());
			int k = sam2[*it];
			model2->addConstr(vars2[x] + vars2[k + mi.size()], GRB_LESS_EQUAL, 1);
		}
	}


	// shadj coexist constraints
	sp.clear();
	for(mit = mi.begin(); mit != mi.end(); mit++)
	{
		int g = gi[mit->first];
		set<int> s = set_intersection(vct[g], sas2);
		set<int>::iterator it1, it2;
		for(it1 = s.begin(); it1 != s.end(); it1++)
		{
			int x = *it1;
			shadj & sax = vsa[x];
			assert(sam2.find(x) != sam2.end());
			int xi = sam2[x];
			for(it2 = s.begin(); it2 != s.end(); it2++)
			{
				int y = *it2;
				shadj & say = vsa[y];
				assert(sam2.find(y) != sam2.end());
				int yi = sam2[y];

				if(x >= y) continue;
				bool b = sax.coexist(say);
				if(b == true) continue;

				if(sp.find(pair<int, int>(x, y)) != sp.end()) continue;
				sp.insert(pair<int, int>(x, y));

				model2->addConstr(vars2[mi.size() + xi] + vars2[mi.size() + yi], GRB_LESS_EQUAL, 1);
			}
		}
	}


	// objective function
	GRBLinExpr obj2;
	for(int i = 0; i < sav2.size(); i++) obj2 += vars2[i + mi.size()];

	int nsa = xv.size() - 1;

	// optimize
	model2->getEnv().set(GRB_IntParam_OutputFlag, 0);
	model2->setObjective(obj2, GRB_MAXIMIZE);
	model2->getEnv().set(GRB_IntParam_SolutionLimit, 1);
	model2->getEnv().set(GRB_DoubleParam_Cutoff, nsa + 0.5);
	model2->update();

	model2->optimize();
	int code2 = model2->get(GRB_IntAttr_Status);
	int ans2 = model2->get(GRB_DoubleAttr_ObjVal);
	if(code2 == GRB_CUTOFF) ans2 = -1;


	printf("seed length = %3d, %4d affected shadjs, %4d affected genes, %4d variables, opt = %3d, nsa = %3d, %s\n",
			(int)xv.size(), (int)sav2.size(), (int)mi.size(), (int)vars2.size(), ans2, nsa, ans2 <= nsa ? "TRUE" : "FALSE");

	printf("shared adjacencies = ");
	for(int i = 0; i < v0.size(); i++) printf("%5d,", v0[i]);
	printf("\n");

	printf("segment1: (");
	for(int i = 0; i < xv.size(); i++) printf("%4d[%4d], ", gi[xv[i]], xv[i]->x);
	printf(")\n");

	printf("segment2: (");
	for(int i = 0; i < yv.size(); i++) printf("%4d[%4d], ", gi[yv[i]], yv[i]->x);
	printf(")\n");

	for(int i = 0; i < sav2.size(); i++)
	{
		bool f = true;
		if(ans2 >= 0 && vars2[i + mi.size()].get(GRB_DoubleAttr_X) <= 0) f = false;
		if(f == false) continue;
		string s = vsa[sav2[i]].print_string(gi);
		printf("sh-adjs(%4d): %s [%s]\n", sav2[i], s.c_str(), f ? "true" : "false");
	}
	printf("\n");

	delete model2;
	
	if(ans2 >= nsa + 1) return false;

	for(int k = 0; k < xv.size(); k++) add_fixed_pair(xv[k], yv[k]);
	return true;

	// collect inference constraints
	/*
	vector<int> v;
	for(int k = 0; k < sav2.size(); k++)
	{
		if(vars2[k + mi.size()].get(GRB_DoubleAttr_X) <= 0) continue;
		assert(vars2[k + mi.size()].get(GRB_DoubleAttr_X) >= 1);
		v.push_back(k);
	}

	printf("INFERENCE TESTING\n");
	for(int k = 0; k < v.size(); k++)
	{
		GRBConstr c = model2->addConstr(vars2[k + mi.size()], GRB_EQUAL, 0);
		model2->update();

		model2->optimize();
		code = model2->get(GRB_IntAttr_Status);
		ans2 = model2->get(GRB_DoubleAttr_ObjVal);
		if(code == GRB_CUTOFF) ans2 = -1;

		if(ans2 <= nsa) printf("INFERENCE CONSTRAINTS FOUND\n");

		model2->remove(c);
		model2->update();
	}
	*/
}

int psolver::check(const set<int> & saopt)
{
	for(int i = 0; i < cds.size(); i++)
	{
		vector<int> v = locate_shared_adjacencies(cds[i]);
		int n = 0;
		for(int k = 0; k < v.size(); k++)
		{
			if(saopt.find(v[k]) != saopt.end()) n++;
		}
		assert(n == 0 || n == v.size());
		if(n == 0) continue;

		candidate & cd = cds[i];

		for(int k = 0; k < cd.xv.size(); k++)
		{
			if(is_fixed_gene(cd.xv[k])) continue;
			if(is_fixed_gene(cd.yv[k])) continue;
		}

		verifier vr(cd.xv, cd.yv, cd.d, this, env, conf->max_combination);
		vr.verify();

		char s[1024];
		sprintf(s, "tex/segment%d.tex", i);
		vr.draw(s);
	}

	return 0;
}

bool psolver::remove_redundant_genes()
{
	bool b = false;
	for(int f = 1; f < gf1.size(); f++)
	{
		if(vf[f] == 0) continue;
		for(int k = 0; k < gf1[f].size(); k++)
		{
			gene * g = gf1[f][k];
			if(is_fixed_gene(g) == true) continue;
			if(vct[gi[g]].size() >= 1) continue;
			operation * op = new deletion(g, g, true);
			gm1->operate(op);
			delete op;
			b = true;
		}
		for(int k = 0; k < gf2[f].size(); k++)
		{
			gene * g = gf2[f][k];
			if(is_fixed_gene(g) == true) continue;
			if(vct[gi[g]].size() >= 1) continue;
			operation * op = new deletion(g, g, true);
			gm2->operate(op);
			delete op;
			b = true;
		}
	}
	return b;
}
