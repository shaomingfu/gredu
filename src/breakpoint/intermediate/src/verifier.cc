#include "verifier.h"

#include <algorithm>
#include "draw.h"

verifier::verifier(const vector<gene*> & x, const vector<gene*> & y, bool dd, pbase * pp, GRBEnv * e, int max)
	: xv(x), yv(y), d(dd), p(pp), env(e), max_combination(max)
{
	v0 = p->locate_shared_adjacencies(xv, yv);
	feasible = true;
	ans1 = 0;
	ans2 = 9999;
}

verifier::~verifier()
{}

bool verifier::verify()
{
	build_affected_shared_adjacencies();
	if(feasible == false) return false;

	build_affected_genes();
	if(feasible == false) return false;

	build_innocent_shared_adjacencies();
	if(feasible == false) return false;
	
	optimize_innocent_shared_adjacencies();
	if(feasible == false) return false;

	optimize_affected_shared_adjacencies();
	if(feasible == false) return false;

	//print_optimal_solutions();
	
	if(ans2 >= (int)xv.size()) return false;

	for(int k = 0; k < xv.size(); k++)
	{
		p->add_fixed_pair(xv[k], yv[k]);
	}

	return true;
}

int verifier::build_affected_shared_adjacencies()
{
	int n = 0;
	for(int i = 0; i < xv.size(); i++)
	{
		int x = p->vct[p->gi[xv[i]]].size();
		if(n >= max_combination) feasible = false;
		if(feasible == false) return 0;
		n += x;
	}
	for(int i = 0; i < yv.size(); i++)
	{
		int x = p->vct[p->gi[yv[i]]].size();
		if(n >= max_combination) feasible = false;
		if(feasible == false) return 0;
		n += x;
	}

	if(n >= 2 * max_combination) feasible = false;
	if(feasible == false) return 0;

	set<int> spx = p->build_span_intersection(xv);
	set<int> spy = p->build_span_intersection(yv);

	set<int> ctx = p->build_contact_union(xv);
	set<int> cty = p->build_contact_union(yv);

	set<int> ssa;
	ssa = set_union(ssa, spx);
	ssa = set_union(ssa, spy);
	ssa = set_union(ssa, ctx);
	ssa = set_union(ssa, cty);

	sam2.clear();
	sav2.clear();
	// construct affected shared adjacencies
	set<int>::iterator it;
	for(it = ssa.begin(); it != ssa.end(); it++)
	{
		shadj & sa = p->vsa[*it];
		if(p->innocent(sa, xv, yv) == true) continue;
		sam2.insert(PI(*it, sav2.size()));
		sav2.push_back(*it);
	}

	return 0;
}

int verifier::build_affected_genes()
{
	// construct affected genes
	mi.clear();
	for(int i = 0; i < xv.size(); i++)
	{
		mi.insert(pair<gene*, int>(xv[i], mi.size()));
		mi.insert(pair<gene*, int>(yv[i], mi.size()));
	}

	for(int i = 0; i < sav2.size(); i++)
	{
		shadj & sa = p->vsa[sav2[i]];
		vector<gene*> v = sa.get_genes();
		for(int k = 0; k < v.size(); k++)
		{
			mi.insert(pair<gene*, int>(v[k], mi.size()));
		}
	}

	if(mi.size() >= max_combination) feasible = false;

	return 0;
}

int verifier::build_innocent_shared_adjacencies()
{
	set<int> ssa;
	map<gene*, int>::iterator mit;
	for(mit = mi.begin(); mit != mi.end(); mit++)
	{
		set<int> & s = p->vct[p->gi[mit->first]];
		ssa = set_union(ssa, s);
	}

	// construct innocent shared adjacencies
	sam1.clear();
	sav1.clear();
	set<int>::iterator it;
	for(it = ssa.begin(); it != ssa.end(); it++)
	{
		shadj & sa = p->vsa[*it];

		//if(p->innocent(sa, xv, yv) == false) continue;

		if(mi.find(sa.x1) == mi.end()) continue;
		if(mi.find(sa.x2) == mi.end()) continue;
		if(mi.find(sa.y1) == mi.end()) continue;
		if(mi.find(sa.y1) == mi.end()) continue;

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

	return 0;
}

int verifier::optimize_innocent_shared_adjacencies()
{
	if(mi.size() + sav1.size() >= max_combination) feasible = false;
	if(feasible == false) return 0;

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
		shadj sa = p->vsa[sav1[i]];
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
		int g = p->gi[mit->first];
		int x = mit->second;
		set<int> s = set_intersection(p->vsp[g], sas1);

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
		int g = p->gi[mit->first];
		set<int> s = set_intersection(p->vct[g], sas1);
		set<int>::iterator it1, it2;
		for(it1 = s.begin(); it1 != s.end(); it1++)
		{
			int x = *it1;
			shadj & sax = p->vsa[x];
			assert(sam1.find(x) != sam1.end());
			int xi = sam1[x];
			for(it2 = s.begin(); it2 != s.end(); it2++)
			{
				int y = *it2;
				shadj & say = p->vsa[y];
				assert(sam1.find(y) != sam1.end());
				int yi = sam1[y];

				if(x >= y) continue;
				bool b = sax.coexist(say);
				if(b == true) continue;

				if(sp.find(pair<int, int>(x, y)) != sp.end()) continue;
				sp.insert(pair<int, int>(x, y));

				//printf("xi = %5d, yi = %5d\n", sav1[xi], sav1[yi]);
				model1->addConstr(vars1[mi.size() + xi] + vars1[mi.size() + yi], GRB_LESS_EQUAL, 1);
			}
		}
	}

	/* testing
	for(int i = 0; i < v0.size(); i++)
	{
		string s = p->vsa[v0[i]].print_string(p->gi);
		printf("given shadj = %5d: %s\n", v0[i], s.c_str());
	}

	for(int i = 0; i < xv.size(); i++)
	{
		printf("xv%5d = %6d (%5d) || yv%5d = %6d (%5d)\n", i, xv[i]->x, p->gi[xv[i]], i, yv[i]->x, p->gi[yv[i]]);
	}
	*/

	// objective function
	GRBLinExpr obj1;
	for(int i = 0; i < sav1.size(); i++) obj1 += vars1[i + mi.size()];

	// optimize
	model1->getEnv().set(GRB_IntParam_OutputFlag, 0);
	model1->setObjective(obj1, GRB_MAXIMIZE);
	model1->update();

	model1->optimize();

	ans1 = model1->get(GRB_DoubleAttr_ObjVal);

	opt1.clear();
	for(int i = 0; i < sav1.size(); i++)
	{
		if(ans1 >= 0 && vars1[i + mi.size()].get(GRB_DoubleAttr_X) <= 0) continue;
		opt1.push_back(sav1[i]);
	}

	delete model1;

	//printf("given %3d shared adjacencies, optimal value is %3d\n", (int)v0.size(), ans1);

	return 0;
}

int verifier::optimize_affected_shared_adjacencies()
{
	if(mi.size() + sav2.size() >= max_combination) feasible = false;
	if(feasible == false) return 0;

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

	// IMPORTANT
	// for each optimal shared adjacencies in opt1
	for(int i = 0; i < opt1.size(); i++)
	{
		GRBVar var = model2->addVar(0, 1, -1, GRB_BINARY);
		vars2.push_back(var);
	}

	// for each gene indicate whether it is used
	for(int i = 0; i < mi.size(); i++)
	{
		GRBVar var = model2->addVar(0, 1, 0, GRB_BINARY);
		vars2.push_back(var);
	}


	model2->update();

	// Constraints
	// shadj --> boundary gene
	for(int i = 0; i < sav2.size(); i++)
	{
		shadj & sa = p->vsa[sav2[i]];
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
	map<gene*, int>::iterator mit;
	set<int> sas2(sav2.begin(), sav2.end());
	for(mit = mi.begin(); mit != mi.end(); mit++)
	{
		int g = p->gi[mit->first];
		int x = mit->second;
		set<int> s = set_intersection(p->vsp[g], sas2);

		for(set<int>::iterator it = s.begin(); it != s.end(); it++)
		{
			assert(sam2.find(*it) != sam2.end());
			int k = sam2[*it];
			model2->addConstr(vars2[x] + vars2[k + mi.size()], GRB_LESS_EQUAL, 1);
		}
	}

	// shadj coexist constraints
	set<PI> sp;
	for(mit = mi.begin(); mit != mi.end(); mit++)
	{
		int g = p->gi[mit->first];
		set<int> s = set_intersection(p->vct[g], sas2);
		set<int>::iterator it1, it2;
		for(it1 = s.begin(); it1 != s.end(); it1++)
		{
			int x = *it1;
			shadj & sax = p->vsa[x];
			assert(sam2.find(x) != sam2.end());
			int xi = sam2[x];
			for(it2 = s.begin(); it2 != s.end(); it2++)
			{
				int y = *it2;
				shadj & say = p->vsa[y];
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

	// IMPORTANT
	// infer used genes
	for(int i = 0; i < sav2.size(); i++)
	{
		shadj & sa = p->vsa[sav2[i]];
		vector<gene*> v = sa.get_genes();
		int n = mi.size() + sav2.size() + opt1.size();
		for(int k = 0; k < v.size(); k++)
		{
			int x = mi[v[k]];
			model2->addConstr(vars2[n + x], GRB_GREATER_EQUAL, vars2[mi.size() + i]);
		}
	}

	// infer opt1
	for(int i = 0; i < opt1.size(); i++)
	{
		shadj & sa = p->vsa[opt1[i]];

		GRBLinExpr expr;
		int n = mi.size() + sav2.size() + opt1.size();
		expr += vars2[n + mi[sa.x1]];
		expr += vars2[n + mi[sa.x2]];
		expr += vars2[n + mi[sa.y1]];
		expr += vars2[n + mi[sa.y2]];

		int x = mi.size() + sav2.size() + i;
		model2->addConstr(vars2[x], GRB_GREATER_EQUAL, expr - 3);
	}

	// objective function
	GRBLinExpr obj2;
	for(int i = 0; i < sav2.size(); i++) obj2 += vars2[i + mi.size()];
	for(int i = 0; i < opt1.size(); i++) obj2 -= vars2[i + mi.size() + sav2.size()];

	int nsa = xv.size() - 1;

	// optimize
	model2->getEnv().set(GRB_IntParam_OutputFlag, 0);
	model2->setObjective(obj2, GRB_MAXIMIZE);
	model2->getEnv().set(GRB_IntParam_SolutionLimit, 1);
	model2->getEnv().set(GRB_DoubleParam_Cutoff, nsa + 0.5);
	model2->update();


	model2->optimize();
	int code2 = model2->get(GRB_IntAttr_Status);
	ans2 = model2->get(GRB_DoubleAttr_ObjVal);
	if(code2 == GRB_CUTOFF) ans2 = -1;
	
	opt2.clear();
	for(int i = 0; i < sav2.size(); i++)
	{
		if(ans2 >= 0 && vars2[i + mi.size()].get(GRB_DoubleAttr_X) <= 0) continue;
		opt2.push_back(sav2[i]);
	}

	delete model2;
	return 0;
}

int verifier::print_optimal_solutions()
{
	printf("length = %3d, %4d affected shadjs, %4d affected genes, %4d innocent shadjs, ans1 = %3d, ans2 = %3d: %s\n",
			(int)xv.size(), (int)sav2.size(), (int)mi.size(), (int)sav1.size(), ans1, ans2, ans2 < (int)xv.size() ? "TRUE" : "FALSE");

	printf("shared adjacencies = ");
	for(int i = 0; i < v0.size(); i++) printf("%5d,", v0[i]);
	printf("\n");

	printf("segment1: (");
	for(int i = 0; i < xv.size(); i++) printf("%4d[%4d], ", p->gi[xv[i]], xv[i]->x);
	printf(")\n");

	printf("segment2: (");
	for(int i = 0; i < yv.size(); i++) printf("%4d[%4d], ", p->gi[yv[i]], yv[i]->x);
	printf(")\n");

	printf("optimal solution 1:\n");
	set<int>::iterator it;
	for(int i = 0; i < opt1.size(); i++)
	{
		int x = opt1[i];	
		string s = p->vsa[x].print_string(p->gi);
		printf("sh-adjs(%4d): %s\n", x, s.c_str());
	}

	printf("optimal solution 2:\n");
	for(int i = 0; i < opt2.size(); i++)
	{
		int x = opt2[i];	
		string s = p->vsa[x].print_string(p->gi);
		printf("sh-adjs(%4d): %s\n", x, s.c_str());
	}

	printf("\n");

	return 0;
}

int verifier::draw(const string & file)
{
	ofstream fout(file.c_str());
	
	draw_header(fout);

	fout<<"\\def\\glen{1cm}\n";
	fout<<"\\def\\ylen{2cm}\n";

	set<int> sx;
	set<int> sy;
	for(int i = 0; i < opt1.size(); i++)
	{
		shadj & sa = p->vsa[opt1[i]];
		if(sx.find(p->gi[sa.x1]) == sx.end()) sx.insert(p->gi[sa.x1]);
		if(sx.find(p->gi[sa.x2]) == sx.end()) sx.insert(p->gi[sa.x2]);
		if(sy.find(p->gi[sa.y1]) == sy.end()) sy.insert(p->gi[sa.y1]);
		if(sy.find(p->gi[sa.y2]) == sy.end()) sy.insert(p->gi[sa.y2]);
	}

	for(int i = 0; i < opt2.size(); i++)
	{
		shadj & sa = p->vsa[opt2[i]];
		if(sx.find(p->gi[sa.x1]) == sx.end()) sx.insert(p->gi[sa.x1]);
		if(sx.find(p->gi[sa.x2]) == sx.end()) sx.insert(p->gi[sa.x2]);
		if(sy.find(p->gi[sa.y1]) == sy.end()) sy.insert(p->gi[sa.y1]);
		if(sy.find(p->gi[sa.y2]) == sy.end()) sy.insert(p->gi[sa.y2]);
	}


	vector<int> vx(sx.begin(), sx.end());
	vector<int> vy(sy.begin(), sy.end());

	sort(vx.begin(), vx.end());
	sort(vy.begin(), vy.end());
	if(d == false) reverse(vy.begin(), vy.end());

	map<int, int> mx;
	map<int, int> my;
	for(int k = 0; k < vx.size(); k++) mx.insert(PI(vx[k], k));
	for(int k = 0; k < vy.size(); k++) my.insert(PI(vy[k], k));

	int px = -1;
	int py = -1;
	for(int k = 0; k < vx.size(); k++)
	{
		if(p->gi[xv[0]] == vx[k])
		{
			px = k;
			break;
		}
	}
	assert(px != -1);

	for(int k = 0; k < vy.size(); k++)
	{
		if(p->gi[yv[0]] == vy[k])
		{
			py = k;
			break;
		}
	}
	assert(py != -1);
	
	// vertices
	for(int i = 0; i < vx.size(); i++)
	{
		string cx = "colx";
		fout<<"\\node[myrectangle, \\"<<cx.c_str()<<", draw, label = above:"<<vx[i]<<"] ";
		fout<<"(x"<<i<<") at ("<< i - px <<" * \\glen, "<< 0 <<" * \\ylen) {"<< p->ig[vx[i]]->x <<"};\n";
	}

	for(int i = 0; i < vy.size(); i++)
	{
		string cx = "colx";
		fout<<"\\node[myrectangle, \\"<<cx.c_str()<<", draw, label = below:"<<vy[i]<<"] ";
		fout<<"(y"<<i<<") at ("<< i - py <<" * \\glen, "<< -1 <<" * \\ylen) {"<< p->ig[vy[i]]->x <<"};\n";
	}

	// highlight given shared segment
	for(int i = 0; i < xv.size(); i++)
	{
		int x = mx[p->gi[xv[i]]];
		string cx = "colx";
		fout<<"\\node[myrectangle, very thick, \\"<<cx.c_str()<<", draw] ";
		fout<<" at ("<< x - px <<" * \\glen, "<< 0 <<" * \\ylen) {"<< xv[i]->x <<"};\n";
		fout<<" at ("<< x - px <<" * \\glen, "<< 0 <<" * \\ylen) {};\n";
	}

	for(int i = 0; i < yv.size(); i++)
	{
		int x = my[p->gi[yv[i]]];
		string cx = "colx";
		fout<<"\\node[myrectangle, very thick, \\"<<cx.c_str()<<", draw] ";
		fout<<" at ("<< x - py <<" * \\glen, "<< -1 <<" * \\ylen) {"<< yv[i]->x <<"};\n";
	}

	// edges in opt1
	for(int i = 0; i < opt1.size(); i++)
	{
		shadj & sa = p->vsa[opt1[i]];
		int x1 = mx[p->gi[sa.x1]];
		int x2 = mx[p->gi[sa.x2]];
		int y1 = my[p->gi[sa.y1]];
		int y2 = my[p->gi[sa.y2]];
		
		int xx1 = x1 < x2 ? x1 : x2;
		int xx2 = x1 > x2 ? x1 : x2;
		int yy1 = y1 < y2 ? y1 : y2;
		int yy2 = y1 > y2 ? y1 : y2;

		fout<<"\\draw[line width = 0.06cm, \\cola, bend right = 30] (x"<<xx2<<") to (x"<<xx1<<");\n";
		fout<<"\\draw[line width = 0.06cm, \\cola, bend right = 30] (y"<<yy1<<") to (y"<<yy2<<");\n";

		fout<<"\\draw[line width = 0.06cm, \\cola, <-> ] (x"<<x1<<") -- (y"<<y1<<");\n";
		fout<<"\\draw[line width = 0.06cm, \\cola, <-> ] (x"<<x2<<") -- (y"<<y2<<");\n";
	}

	// edges in opt2
	for(int i = 0; i < opt2.size(); i++)
	{
		shadj & sa = p->vsa[opt2[i]];
		int x1 = mx[p->gi[sa.x1]];
		int x2 = mx[p->gi[sa.x2]];
		int y1 = my[p->gi[sa.y1]];
		int y2 = my[p->gi[sa.y2]];
		
		int xx1 = x1 < x2 ? x1 : x2;
		int xx2 = x1 > x2 ? x1 : x2;
		int yy1 = y1 < y2 ? y1 : y2;
		int yy2 = y1 > y2 ? y1 : y2;

		fout<<"\\draw[thick, \\colb, bend right = 30] (x"<<xx2<<") to (x"<<xx1<<");\n";
		fout<<"\\draw[thick, \\colb, bend right = 30] (y"<<yy1<<") to (y"<<yy2<<");\n";

		fout<<"\\draw[thick, \\colb, <-> ] (x"<<x1<<") -- (y"<<y1<<");\n";
		fout<<"\\draw[thick, \\colb, <-> ] (x"<<x2<<") -- (y"<<y2<<");\n";
	}

	// string
	char s[10240];
	sprintf(s, "ans1 = %d, ans2 = %d\n", ans1, ans2);
	fout<<"\\node[myrectangle, \\colc] at (1.5 * \\glen, -1.8 * \\ylen) {"<< s <<"};\n";

	draw_footer(fout);
	fout.close();
	return 0;
}


