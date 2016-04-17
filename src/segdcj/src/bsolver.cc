#include "bsolver.h"
#include "draw.h"
#include "lpsolver2.h"

#include <iostream>
#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <cmath>

bsolver::bsolver(config * _conf, genome * _gm1, genome * _gm2)
	: asolver(_conf, _gm1, _gm2)
{
	return;
	build_possibly_duplicated_genes();
	build_matching_graph();
	build_initial_matching();

	build_shared_segments();
	make_shared_segments();
	prepare_shared_segments();

	int n = verify_shared_segments();
	printf("identify %4d shared segments\n", n);
}

bsolver::~bsolver()
{}

int bsolver::build_shared_segments()
{
	vss.clear();
	edge_iterator ei1, ei2;
	for(tie(ei1, ei2) = edges(ar); ei1 != ei2; ei1++)
	{
		int x = source(*ei1, ar);
		int y = target(*ei1, ar);
		add_shared_segment(x, y);
	}
	return 0;
}

int bsolver::add_shared_segment(int x, int y)
{
	for(int i = 0; i < vss.size(); i++)
	{
		if(vss[i].contain(x / 2, y / 2)) return 0;
	}
		
	vector<int> va;
	vector<int> vb;
	set<int> sx;
	set<int> sy;
	extend_forward(x, y, va, vb, sx, sy);

	vector<int> xv;
	vector<int> yv;
	for(int i = va.size() - 1; i >= 0; i--)
	{
		xv.push_back(va[i]);
		yv.push_back(vb[i]);
	}
	xv.push_back(x / 2);
	yv.push_back(y / 2);

	assert(sx.find(x / 2) == sx.end());
	assert(sy.find(y / 2) == sy.end());
	sx.insert(x / 2);
	sy.insert(y / 2);

	extend_backward(x, y, xv, yv, sx, sy);

	// TODO
	if(xv.size() <= 1) return 0;

	for(int i = 0; i < xv.size(); i++)
	{
		int dx = out_degree(xv[i], pr);
		int dy = out_degree(yv[i], pr);
		if(dx == 1 && dy == 1) return 0;
	}

	int dx = 0;
	int dy = 0;
	for(int i = 0; i < xv.size(); i++)
	{
		dx += out_degree(xv[i], pr);
		dy += out_degree(yv[i], pr);
	}

	assert(dx >= xv.size());
	assert(dy >= yv.size());

	if(dx > xv.size() && dy > yv.size()) return 0;

	// TODO
	if(dx > 3 * dy || dy > 3 * dx) return 0;

	shseg ss(xv, yv);
	vss.push_back(ss);

	return 0;
}

int bsolver::extend_forward(int x, int y, vector<int> & xv, vector<int> & yv, set<int> & sx, set<int> & sy)
{
	assert(edge(x, y, ar).second);
	int s = adj[x];
	int t = adj[y];
	if(s == -1) return 0;
	assert(t != -1);
	if(edge(s, t, ar).second == false) return 0;
	if(sx.find(s / 2) != sx.end()) return 0;
	if(sy.find(t / 2) != sy.end()) return 0;
	xv.push_back(s / 2);
	yv.push_back(t / 2);
	sx.insert(s / 2);
	sy.insert(t / 2);
	extend_backward(s, t, xv, yv, sx, sy);
	return 0;
}

int bsolver::extend_backward(int x, int y, vector<int> & xv, vector<int> & yv, set<int> & sx, set<int> & sy)
{
	assert(edge(x, y, ar).second);
	int s = cpl[x];
	int t = cpl[y];
	assert(edge(s, t, ar).second);
	extend_forward(s, t, xv, yv, sx, sy);
	return 0;
}

int bsolver::build_possibly_duplicated_genes()
{
	for(int i = 0; i < num_vertices(pr); i++)
	{
		if(out_degree(i, pr) <= 1) continue;
		gene * g = (i < gi1.size()) ? ig1[i] : ig2[i - gi1.size()];
		if(g->x == 0) continue;

		out_edge_iterator ei1, ei2;
		for(tie(ei1, ei2) = out_edges(i, pr); ei1 != ei2; ei1++)
		{
			int t = target(*ei1, pr);
			/*
			printf("checking %6d (degree = %5d), %6d. sdup has %6d elements:", i, out_degree(i, pr), t, (int)sdup.size());
			for(set<int>::iterator it = sdup.begin(); it != sdup.end(); it++) printf("%5d ", *it);
			printf("\n");
			*/
			if(sdup.find(t) == sdup.end()) sdup.insert(t);
		}
	}
	return 0;
}

int bsolver::build_matching_graph()
{
	mr.clear();
	int n = num_vertices(ar);
	for(int i = 0; i < 2 * n; i++) add_vertex(mr);

	// add edges for the extremity pairs
	edge_iterator ei1, ei2;
	for(tie(ei1, ei2) = edges(ar); ei1 != ei2; ei1++)
	{
		int s = source(*ei1, ar);
		int t = target(*ei1, ar);
		add_edge(s, t, mr);
	}

	// add edges for the possibly duplicated genes
	set<int>::iterator it;
	for(it = sdup.begin(); it != sdup.end(); it++)
	{
		int g = (*it);
		int x = 2 * g + 0;
		int y = 2 * g + 1;
		add_edge(x, y, mr);
	}

	// add edges for adjacencies
	for(int i = 0; i < adj.size(); i++)
	{
		if(adj[i] == -1) continue;
		if(i < adj[i]) continue;
		add_edge(i + n, adj[i] + n, mr);
	}

	// add bridge edges
	for(int i = 0; i < n; i++) add_edge(i, i + n, mr);

	// index map for the matching graph
	mvim.assign(num_vertices(mr), 0);
	for(int i = 0; i < mvim.size(); i++) mvim[i] = i;
	return 0;
}

int bsolver::build_initial_matching()
{
	int n = num_vertices(mr) / 2;
	mate.assign(2 * n, VNULL);
	for(int i = 0; i < n; i++) mate[i] = i + n;
	for(int i = 0; i < n; i++) mate[i + n] = i;
	assert_mate();
	return 0;
}

int bsolver::make_shared_segments()
{
	for(int i = 0; i < vss.size(); i++)
	{
		make_shared_segment(vss[i]);
	}
	return 0;
}

int bsolver::make_shared_segment(shseg & ss)
{
	ss.sa.clear();
	ss.sb.clear();
	ss.ea.clear();
	ss.eb.clear();

	// build intra-extremities
	assert(ss.xv.size() == ss.yv.size());
	for(int i = 0; i < ss.xv.size(); i++)
	{
		ss.sa.insert(ss.xv[i] * 2 + 0);
		ss.sa.insert(ss.xv[i] * 2 + 1);
		ss.sa.insert(ss.yv[i] * 2 + 0);
		ss.sa.insert(ss.yv[i] * 2 + 1);
	}

	set<int>::iterator it;
	out_edge_iterator ei1, ei2;
	for(it = ss.sa.begin(); it != ss.sa.end(); it++)
	{
		int s = (*it);
		for(tie(ei1, ei2) = out_edges(s, ar); ei1 != ei2; ei1++)
		{
			int t = target(*ei1, ar);
			assert(s != t);
			if(ss.sa.find(t) != ss.sa.end()) continue;
			ss.sa.insert(t);
		}
	}

	// build extra extremities
	for(it = ss.sa.begin(); it != ss.sa.end(); it++)
	{
		int s = (*it);
		int t = adj[s];
		if(t == -1) continue;
		if(ss.sa.find(t) != ss.sa.end()) continue;
		ss.sb.insert(t);
	}

	// compute the boundary of the central segment
	ss.xb.first = ss.xb.second = -1;
	for(int i = 0; i < ss.xv.size(); i++)
	{
		int a = ss.xv[i] * 2 + 0;
		int b = adj[a];
		if(ss.sa.find(b) == ss.sa.end())
		{
			assert(ss.sb.find(b) != ss.sb.end());
			if(ss.xb.first == -1) ss.xb.first = b;
			else if(ss.xb.second == -1) ss.xb.second = b;
			else assert(false);
		}
		a = ss.xv[i] * 2 + 1;
		b = adj[a];
		if(ss.sa.find(b) == ss.sa.end())
		{
			assert(ss.sb.find(b) != ss.sb.end());
			if(ss.xb.first == -1) ss.xb.first = b;
			else if(ss.xb.second == -1) ss.xb.second = b;
			else assert(false);
		}
	}

	ss.yb.first = ss.yb.second = -1;
	for(int i = 0; i < ss.yv.size(); i++)
	{
		int a = ss.yv[i] * 2 + 0;
		int b = adj[a];
		if(ss.sa.find(b) == ss.sa.end())
		{
			assert(ss.sb.find(b) != ss.sb.end());
			if(ss.yb.first == -1) ss.yb.first = b;
			else if(ss.yb.second == -1) ss.yb.second = b;
			else assert(false);
		}
		a = ss.yv[i] * 2 + 1;
		b = adj[a];
		if(ss.sa.find(b) == ss.sa.end())
		{
			assert(ss.sb.find(b) != ss.sb.end());
			if(ss.yb.first == -1) ss.yb.first = b;
			else if(ss.yb.second == -1) ss.yb.second = b;
			else assert(false);
		}
	}

	// build intra-edges
	for(it = ss.sa.begin(); it != ss.sa.end(); it++)
	{
		int s = (*it);
		for(tie(ei1, ei2) = out_edges(s, ar); ei1 != ei2; ei1++)
		{
			int t = target(*ei1, ar);
			assert(s != t);
			if(s > t) continue;
			assert(s < gi1.size() * 2);
			PI p(s, t);
			if(ss.ea.find(p) != ss.ea.end()) continue;
			ss.ea.insert(p);
		}
	}

	// build extra edges
	// remove intra-edges from the matching graph
	set<PI>::iterator ei;
	for(ei = ss.ea.begin(); ei != ss.ea.end(); ei++)
	{
		int x = ei->first;
		int y = ei->second;
		remove_edge(x, y, mr);
	}
	for(set<int>::iterator it = ss.sa.begin(); it != ss.sa.end(); it++)
	{
		int x = (*it);
		int y = x % 2 == 0 ? x + 1 : x - 1;
		assert(ss.sa.find(y) != ss.sa.end());
		if(x > y) continue;
		int g = x / 2;
		if(sdup.find(g) == sdup.end()) continue;
		assert(edge(x, y, mr).second == true);
		remove_edge(x, y, mr);
	}

	// check possible extra edges
	ss.eb.clear();
	set<int>::iterator it1;
	set<int>::iterator it2;
	for(it1 = ss.sb.begin(); it1 != ss.sb.end(); it1++)
	{
		int x = (*it1);
		for(it2 = ss.sb.begin(); it2 != ss.sb.end(); it2++)
		{
			int y = (*it2);
			if(x <= y) continue;
			//if(x >= gi1.size() * 2) continue;
			//if(y < gi1.size() * 2) continue;
			bool b = check_alternating_path(x, y);
			if(b == false) continue;
			PI p(x, y);
			assert(ss.eb.find(p) == ss.eb.end());
			ss.eb.insert(p);
		}
	}
	
	// recover those edges for matching graph
	for(ei = ss.ea.begin(); ei != ss.ea.end(); ei++)
	{
		int x = ei->first;
		int y = ei->second;
		add_edge(x, y, mr);
	}
	for(set<int>::iterator it = ss.sa.begin(); it != ss.sa.end(); it++)
	{
		int x = (*it);
		int y = x % 2 == 0 ? x + 1 : x - 1;
		assert(ss.sa.find(y) != ss.sa.end());
		if(x > y) continue;
		int g = x / 2;
		if(sdup.find(g) == sdup.end()) continue;
		assert(edge(x, y, mr).second == false);
		add_edge(x, y, mr);
	}

	return 0;
}

bool bsolver::check_alternating_path(int i, int j)
{
	int n = num_vertices(ar);
	clear_vertex(i + n, mr);
	clear_vertex(j + n, mr);
	mate[i] = VNULL;
	mate[j] = VNULL;
	mate[i + n] = VNULL;
	mate[j + n] = VNULL;

	//const_vertex_index_map vim = get(vertex_index, mr);
	path_finder finder(mr, &mate[0], &mvim[0]);
	bool b = finder.augment_matching();

	mate[i] = i + n;
	mate[j] = j + n;
	mate[i + n] = i;
	mate[j + n] = j;
	add_edge(i, i + n, mr);
	add_edge(j, j + n, mr);
	add_edge(i + n, adj[i] + n, mr);
	add_edge(j + n, adj[j] + n, mr);

	assert_mate();
	return b;
}

int bsolver::prepare_shared_segments()
{
	for(int i = 0; i < vss.size(); i++)
	{
		prepare_shared_segment(vss[i]);
	}
	return 0;
}

int bsolver::prepare_shared_segment(shseg & ss)
{
	// build vertex map e2i
	ss.ei.clear();
	ss.ie.clear();
	for(set<int>::iterator it = ss.sa.begin(); it != ss.sa.end(); it++)
	{
		int e = (*it);
		int i = ss.ei.size();
		ss.ei.insert(PI(e, i));
		ss.ie.insert(PI(i, e));
	}
	for(set<int>::iterator it = ss.sb.begin(); it != ss.sb.end(); it++)
	{
		int e = (*it);
		int i = ss.ei.size();
		ss.ei.insert(PI(e, i));
		ss.ie.insert(PI(i, e));
	}

	// build adjacency graph
	ss.ar.clear();
	for(int i = 0; i < ss.ei.size(); i++) add_vertex(ss.ar);
	for(set<PI>::iterator it = ss.ea.begin(); it != ss.ea.end(); it++)
	{
		// remove those edges in the central segment
		if(ss.contain(it->first / 2, it->second / 2)) continue;
		int x = ss.ei[it->first];
		int y = ss.ei[it->second];
		add_edge(x, y, ss.ar);
	}
	for(set<PI>::iterator it = ss.eb.begin(); it != ss.eb.end(); it++)
	{
		int x = ss.ei[it->first];
		int y = ss.ei[it->second];
		add_edge(x, y, ss.ar);
	}
	
	// build list of genes
	ss.ge.clear();
	ss.gf.clear();
	for(set<int>::iterator it = ss.sa.begin(); it != ss.sa.end(); it++)
	{
		int x = (*it);
		int y = (x % 2 == 0) ? x + 1 : x - 1;
		assert(ss.sa.find(y) != ss.sa.end());
		if(x % 2 == 1) continue;
		ss.ge.push_back(PI(ss.ei[x], ss.ei[y]));
		int t = x / 2;
		if(t < ig1.size())
		{
			gene * g = ig1[t];
			int f = (int)fabs(g->x);
			ss.gf.push_back(f);
		}
		else
		{
			gene * g = ig2[t - ig1.size()];
			int f = (int)fabs(g->x);
			ss.gf.push_back(0 - f);
		}
	}

	/*
	for(set<int>::iterator it = ss.sb.begin(); it != ss.sb.end(); it++)
	{
		int x = (*it);
		int y = (x % 2 == 0) ? x + 1 : x - 1;
		if(ss.sb.find(y) == ss.sb.end()) continue;
		if(x % 2 == 1) continue;
		ss.ge.push_back(PI(ss.ei[x], ss.ei[y]));
		int t = x / 2;
		if(t < ig1.size())
		{
			gene * g = ig1[t];
			int f = (int)fabs(g->x);
			ss.gf.push_back(f);
		}
		else
		{
			gene * g = ig2[t - ig1.size()];
			int f = (int)fabs(g->x);
			ss.gf.push_back(0 - f);
		}
	}
	*/

	// build adjacent extremties
	ss.adj.assign(ss.ei.size(), -1);
	for(map<int, int>::iterator it = ss.ei.begin(); it != ss.ei.end(); it++)
	{
		int e1 = it->first;
		int i1 = it->second;
		int e2 = adj[e1];
		assert(ss.ei.find(e2) != ss.ei.end());
		int i2 = ss.ei[e2];
		ss.adj[i1] = i2;
	}
	for(int i = 0; i < ss.adj.size(); i++)
	{
		int j = ss.adj[i];
		assert(j < ss.adj.size());
		assert(ss.adj[j] == i);
	}

	// build duplicons (TODO, now assume all are single-gene duplications)
	map< int, vector<int> > m;
	for(int i = 0; i < ss.gf.size(); i++)
	{
		int f = ss.gf[i];
		if(m.find(f) == m.end())
		{
			vector<int> v;
			v.push_back(i);
			m.insert(pair< int, vector<int> >(f, v));
		}
		else
		{
			m[f].push_back(i);
		}
	}

	ss.vdup.clear();
	for(map<int, vector<int> >::iterator it = m.begin(); it != m.end(); it++)
	{
		vector<int> & v = it->second;
		if(v.size() <= 1) continue;
		for(int k = 0; k < v.size(); k++)
		{
			vector<int> dup;
			dup.push_back(v[k]);
			ss.vdup.push_back(dup);
		}
	}
	return 0;
}

int bsolver::verify_shared_segments()
{
	sort(vss.begin(), vss.end());
	int cnt = 0;
	for(int i = 0; i < vss.size(); i++)
	{
		bool b = verify_shared_segment(vss[i]);
		if(b == true) fix_shared_segment(vss[i]);
		if(b == true) cnt++;
	}
	resolve_removed_genes();
	return cnt;
}

bool bsolver::verify_shared_segment(shseg & ss)
{
	lpsolver2 lp(conf, ss.ar, ss.ge, ss.gf, ss.adj, ss.vdup, 0);
	lp.solve();
	int upper = lp.get_num_cycles();
	int lower = ss.calc_lower_bound();

	draw_adjacency_graph("graph.tex", ss);
	ss.draw("shseg.tex", gi1.size() * 2);
	printf("verify = [%c]: upper = %3d, lower = %3d. segment has %2d genes, %3d + %3d extremities, %4d + %4d edges.\n", 
			lower >= upper ? 'T' : 'F', upper, lower, (int)ss.xv.size(), (int)ss.sa.size(), (int)ss.sb.size(), (int)ss.ea.size(), (int)ss.eb.size());

	if(lower >= upper) return true;
	else return false;
}

int bsolver::fix_shared_segment(const shseg & ss)
{
	for(int i = 0; i < ss.xv.size(); i++)
	{
		int x = ss.xv[i];
		int y = ss.yv[i];
		fix_gene_pair(x, y);
		fix_extremity_pair(x * 2 + 0, y * 2 + 0);
		fix_extremity_pair(x * 2 + 1, y * 2 + 1);
	}
	return 0;
}

int bsolver::statistic()
{
	printf("total %5d edges, %5d shared segments\n", (int)num_edges(pr), (int)vss.size());
	for(int i = 0; i < vss.size(); i++) vss[i].statistic();
	return 0;
}

int bsolver::print_shared_segments()
{
	char file1[1024];
	char file2[1024];
	for(int i = 0; i < vss.size(); i++)
	{
		vss[i].print();
		sprintf(file1, "graph%d.tex", i);	
		sprintf(file2, "shseg%d.tex", i);	
		draw_adjacency_graph(file1, vss[0]);
		vss[i].draw(file2, gi1.size() * 2);
	}
	return 0;
}

int bsolver::assert_mate()
{
	int n = num_vertices(ar);
	assert(mate.size() == 2 * n);
	for(int i = 0; i < n; i++) assert(mate[i] == i + n);
	for(int i = 0; i < n; i++) assert(mate[i + n] == i);
	return 0;
}

int bsolver::draw_adjacency_graph(const string & file, shseg & ss)
{
	ofstream fout(file.c_str());
	
	draw_header(fout);

	fout<<"\\def\\glen{1cm}\n";
	fout<<"\\def\\ylen{6cm}\n";

	// vertices
	char sx[1024];
	char sy[1024];
	for(int i = 0; i < gi1.size(); i++)
	{
		int a = i * 2 + 0;
		int b = i * 2 + 1;
		sprintf(sx, "s%d", a);
		sprintf(sy, "s%d", b);
		int x = -1;
		int y = -1;
		if(ig1[i]->x > 0 || ig1[i]->a == NULL)
		{
			x = i * 2 + 0;
			y = i * 2 + 1;
		}
		else
		{
			x = i * 2 + 1;
			y = i * 2 + 0;
		}
		string cx = "colx";
		string cy = "colx";
		if(ss.sa.find(a) != ss.sa.end()) cx = "colc";
		if(ss.sb.find(a) != ss.sb.end()) cx = "colb";
		if(ss.sa.find(b) != ss.sa.end()) cy = "colc";
		if(ss.sb.find(b) != ss.sb.end()) cy = "colb";

		int la = -1;
		int lb = -1;
		if(ss.ei.find(a) != ss.ei.end()) la = ss.ei[a];
		if(ss.ei.find(b) != ss.ei.end()) lb = ss.ei[b];
		fout<<"\\node[mycircle, \\"<<cx.c_str()<<", draw, label = above:"<<la<<"] ("<<sx<<") at ("<<x<<" * \\glen, "<<0<<" * \\ylen) {"<<a<<"};\n";
		fout<<"\\node[mycircle, \\"<<cy.c_str()<<", draw, label = above:"<<lb<<"] ("<<sy<<") at ("<<y<<" * \\glen, "<<0<<" * \\ylen) {"<<b<<"};\n";
	}

	for(int i = 0; i < gi2.size(); i++)
	{
		int a = i * 2 + 0 + gi1.size() * 2;
		int b = i * 2 + 1 + gi1.size() * 2;
		sprintf(sx, "s%d", a);
		sprintf(sy, "s%d", b);
		int x = -1;
		int y = -1;
		if(ig2[i]->x > 0 || ig2[i]->a == NULL)
		{
			x = i * 2 + 0;
			y = i * 2 + 1;
		}
		else
		{
			x = i * 2 + 1;
			y = i * 2 + 0;
		}

		string cx = "colx";
		string cy = "colx";
		if(ss.sa.find(a) != ss.sa.end()) cx = "colc";
		if(ss.sb.find(a) != ss.sb.end()) cx = "colb";
		if(ss.sa.find(b) != ss.sa.end()) cy = "colc";
		if(ss.sb.find(b) != ss.sb.end()) cy = "colb";

		int la = -1;
		int lb = -1;
		if(ss.ei.find(a) != ss.ei.end()) la = ss.ei[a];
		if(ss.ei.find(b) != ss.ei.end()) lb = ss.ei[b];
		fout<<"\\node[mycircle, \\"<<cx.c_str()<<", draw, label = below:"<<la<<"] ("<<sx<<") at ("<<x<<" * \\glen, "<<-1<<" * \\ylen) {"<<a<<"};\n";
		fout<<"\\node[mycircle, \\"<<cy.c_str()<<", draw, label = below:"<<lb<<"] ("<<sy<<") at ("<<y<<" * \\glen, "<<-1<<" * \\ylen) {"<<b<<"};\n";
	}

	// edges in ss.eb
	for(set<PI>::iterator it = ss.eb.begin(); it != ss.eb.end(); it++)
	{
		int x = it->first < it->second ? it->first : it->second;
		int y = it->first > it->second ? it->first : it->second;
		sprintf(sx, "s%d", x);
		sprintf(sy, "s%d", y);
		if(x < gi1.size() * 2 && y >= gi1.size() * 2)
		{
			fout<<"\\draw[double, \\colb] ("<<sx<<") -- ("<<sy<<");\n";
			//fout<<"\\draw[decorate, decoration={coil,aspect=0}, \\colb] ("<<sx<<") -- ("<<sy<<");\n";
		}
		else if(x < gi1.size() * 2 && y < gi1.size() * 2)
		{
			fout<<"\\draw[double, \\colb, bend left = 30] ("<<sx<<") to ("<<sy<<");\n";
			//fout<<"\\draw[decorate, bend left = 30, decoration={coil,aspect=0}, \\colb] ("<<sx<<") to ("<<sy<<");\n";
		}
		else if(x >= gi1.size() * 2 && y >= gi1.size() * 2)
		{
			fout<<"\\draw[double, \\colb, bend right = 30] ("<<sx<<") to ("<<sy<<");\n";
			//fout<<"\\draw[decorate, bend right = 30, decoration={coil,aspect=0}, \\colb] ("<<sx<<") to ("<<sy<<");\n";
		}
		else assert(false);
	}

	// edges between two genomes
	edge_iterator ei1, ei2;
	for(tie(ei1, ei2) = edges(ar); ei1 != ei2; ei1++)
	{
		int x = source(*ei1, ar);
		int y = target(*ei1, ar);
		assert(x < y);
		string c = "colx";
		if(ss.ea.find(PI(x, y)) != ss.ea.end()) c = "colc";
		sprintf(sx, "s%d", x);
		sprintf(sy, "s%d", y);
		string b = "";
		if(ss.contain(x / 2, y / 2)) b = "thick";
		fout<<"\\draw["<<b.c_str()<<",\\"<<c.c_str()<<"] ("<<sx<<") -- ("<<sy<<");\n";
	}
		
	// edges for genes
	for(int i = 0; i < gi1.size(); i++)
	{
		int x = i * 2 + 0;
		int y = i * 2 + 1;
		sprintf(sx, "s%d", x);
		sprintf(sy, "s%d", y);
		string c = "cold";
		if(sdup.find(i) != sdup.end()) c = "colf";
		fout<<"\\draw[line width = 0.07cm, \\"<<c.c_str()<<", ->] ("<<sx<<") to node[label=above:"<<ig1[i]->x<<"]{} ("<<sy<<");\n";
	}

	for(int i = 0; i < gi2.size(); i++)
	{
		int x = i * 2 + 0 + gi1.size() * 2;
		int y = i * 2 + 1 + gi1.size() * 2;
		sprintf(sx, "s%d", x);
		sprintf(sy, "s%d", y);
		string c = "cold";
		if(sdup.find(i + gi1.size()) != sdup.end()) c = "colf";
		fout<<"\\draw[line width = 0.07cm, \\"<<c.c_str()<<", ->] ("<<sx<<") to node[label=below:"<<ig2[i]->x<<"]{} ("<<sy<<");\n";
	}

	// edges for adjacencies
	for(int i = 0; i < adj.size(); i++)
	{
		int j = adj[i];
		if(j == -1) continue;
		if(i > j) continue;
		sprintf(sx, "s%d", i);
		sprintf(sy, "s%d", j);
		if(i >= gi1.size() * 2) fout<<"\\draw[very thick, \\cola, bend left = 30] ("<<sx<<") to ("<<sy<<");\n";
		else fout<<"\\draw[very thick, \\cola, bend right = 30] ("<<sx<<") to ("<<sy<<");\n";
	}

	draw_footer(fout);
	fout.close();
	return 0;
}


