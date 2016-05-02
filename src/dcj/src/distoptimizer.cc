#include "distoptimizer.h"
#include "presolver.h"
#include "ilp6.h"

#include <boost/graph/connected_components.hpp>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <algorithm>

using namespace std;

distoptimizer::distoptimizer(config * _conf, genome * gm1, genome * gm2)
	: conf(_conf)
{
	build_adjacency_sets(gm1, gm2);
	simplified_cycles = 0;
}

distoptimizer::~distoptimizer()
{
	set<gene*>::iterator it;
	for(it = ghosts.begin(); it != ghosts.end(); it++)
	{
		if((*it) != NULL) delete (*it);
	}
}

int distoptimizer::solve()
{
	printf("build graph ...\n");
	build_graph();
	printf("graph has %7d vertices and %7d edges\n", (int)num_vertices(gr), (int)num_edges(gr));

	while(true)
	{
		printf("simplify graph ...\n");
		presolver ps(gr, gene_list, partners);
		bool b = ps.simplify(2);

		simplify_graph();
		printf("graph has %7d vertices and %7d edges\n", (int)num_vertices(gr), (int)num_edges(gr));

		if(b == false) break;
	}
	
	int num1 = 0;
	int num2 = 0;
	int num3 = 0;
	for(int i = 0; i < num_vertices(gr) / 2; i++)
	{
		if(out_degree(i, gr) == 1) num1++;
		if(out_degree(i, gr) == 2) num2++;
		if(out_degree(i, gr) == 3) num3++;
	}
	/*
	printf("%5d out of %5d has degree 1\n", num1, (int)num_vertices(gr) / 2);
	printf("%5d out of %5d has degree 2\n", num2, (int)num_vertices(gr) / 2);
	printf("%5d out of %5d has degree 3\n", num3, (int)num_vertices(gr) / 2);
	printf("\nsimplified_cycles = %5d\n\n", simplified_cycles);
	*/

	printf("run ILP ...\n");

	ilp_base * lp = new ilp6(gr, gene_list, partners);

	lp->set_timelimit(conf->ilp_timelimit);
	lp->model->getEnv().set(GRB_IntParam_MIPFocus, 1);

	lp->solve();
	build_mapping(lp->pairs);
	delete lp;

	int n = sx.size() * 4;
	vector<int> v(n);
	int c = connected_components(bg, &v[0]);
	printf("simplified cycles = %d\n", simplified_cycles);
	printf("\ntotoal %5d adjacencies, total %10d cycles, DCJ distance = %5d\n", n / 4, c, n / 4 - c);
	//print_breakpoint_graph();

	return 0;
}

int distoptimizer::build_adjacency_sets(genome * gm1, genome * gm2)
{
	assert(gm1->alphabet_size == gm2->alphabet_size);

	vector<int> copy1 = gm1->build_gene_copy();
	vector<int> copy2 = gm2->build_gene_copy();
	assert(copy1.size() == copy2.size());

	sx.clear();
	sy.clear();

	gm1->build_adjacencies(sx);
	gm2->build_adjacencies(sy);

	// add ghost genes and adjacencies
	for(int i = 0; i < copy1.size(); i++)
	{
		//assert(copy1.at(i) == copy2.at(i));
		int n = copy1.at(i) - copy2.at(i);

		if(n == 0) continue;
		else if(n > 0)
		{
			for(int k = 0; k < n; k++)
			{
				gene * g = new gene(i + 1, NULL);
				ghosts.insert(g);
				sy.push_back(adjacency(PG(g, g)));
			}
		}
		else if(n < 0)
		{
			n = 0 - n;
			for(int k = 0; k < n; k++)
			{
				gene * g = new gene(i + 1, NULL);
				ghosts.insert(g);
				sx.push_back(adjacency(PG(g, g)));
			}
		}
	}

	int n1 = 0;
	int n2 = 0;
	for(int i = 0; i < gm1->chrms.size(); i++) if(gm1->chrms.at(i)->type == LINEAR) n1++;
	for(int i = 0; i < gm2->chrms.size(); i++) if(gm2->chrms.at(i)->type == LINEAR) n2++;

	int n = n1 - n2;
	if(n > 0)
	{
		for(int k = 0; k < n; k++)
		{
			gene * g1 = new gene(0, NULL);
			gene * g2 = new gene(0, NULL);
			ghosts.insert(g1);
			ghosts.insert(g2);
			sy.push_back(adjacency(PG(g1, g2)));
		}
	}
	else if(n < 0)
	{
		n = 0 - n;
		for(int k = 0; k < n; k++)
		{
			gene * g1 = new gene(0, NULL);
			gene * g2 = new gene(0, NULL);
			ghosts.insert(g1);
			ghosts.insert(g2);
			sx.push_back(adjacency(PG(g1, g2)));
		}
	}

	assert(sx.size() == sy.size());

	return 0;
}

int distoptimizer::build_mapping(const set<PI> & m)
{
	set<PI>::iterator it;
	for(it = m.begin(); it != m.end(); it++)
	{
		int xx = b2a[it->first];
		int yy = b2a[it->second];
		add_extremity_pair(xx, yy);
	}
	return 0;
}

int distoptimizer::add_gene_pair(gene * gx, gene * gy)
{
	map<gene*, gene*>::iterator itx = x2y.find(gx);
	map<gene*, gene*>::iterator ity = y2x.find(gy);

	if(itx == x2y.end())
	{
		assert(ity == y2x.end());
		x2y.insert(PG(gx, gy));
		y2x.insert(PG(gy, gx));
	}
	else
	{
		assert(ity != y2x.end());
		assert(itx->second == gy);
		assert(ity->second == gx);
	}
	return 0;
}

int distoptimizer::add_extremity_pair(int x, int y)
{
	assert(x != y);
	if(x > y)
	{
		int t = x;
		x = y;
		y = t;
	}

	assert(x < sx.size() * 2);
	assert(y >= sx.size() * 2);

	// build the final breakpoint graph
	pair<edge_descriptor, bool> p = edge(x, y, bg);
	if(p.second == false) add_edge(x, y, bg);

	int x1 = (x % 2 == 0) ? (x - 1) : (x + 1);
	int y1 = (y % 2 == 0) ? (y - 1) : (y + 1);
	if(x1 >= 0 && x1 < sx.size() * 2 && y1 >= sx.size() * 2 && y1 < sx.size() * 2 + sy.size() * 2)
	{
		p = edge(x1, y1, bg);
		if(p.second == false) add_edge(x1, y1, bg);
	}

	gene * gx;
	gene * gy;

	if(x % 2 == 0) gx = sx.at(x / 2).e1.g;
	else gx = sx.at(x / 2).e2.g;

	if(y % 2 == 0) gy = sy.at(y / 2 - sx.size()).e1.g;
	else gy = sy.at(y / 2 - sx.size()).e2.g;

	add_gene_pair(gx, gy);
	return 0;
}

int distoptimizer::check_edge(const extremity & ex, const extremity & ey)
{
	if(ex.weak_compare(ey) == 0) return -1;

	map<gene*, gene*>::iterator itx = x2y.find(ex.g);
	map<gene*, gene*>::iterator ity = y2x.find(ey.g);

	if(itx == x2y.end() && ity == y2x.end())
	{
		return 0;
	}
	else if(itx == x2y.end() && ity != x2y.end())
	{
		assert(ity->second != ex.g);
		return -1;
	}
	else if(itx != x2y.end() && ity == x2y.end())
	{
		assert(itx->second != ey.g);
		return -1;
	}
	else if(itx->second == ey.g)
	{
		assert(ity->second == ex.g);
		return 1;
	}
	else
	{
		assert(itx->second != ey.g);
		assert(ity->second != ex.g);
		return -1;
	}
	assert(1 == 0);
	return -1;
}

int distoptimizer::build_graph()
{
	assert(sx.size() == sy.size());
	int n = sx.size() * 2;

	// build the final breakpoint graph
	for(int i = 0; i < 2 * n; i++) add_vertex(bg);
	for(int i = 0; i < n; i++) add_edge(i * 2, i * 2 + 1, bg);

	for(int i = 0; i < 2 * n; i++) add_vertex(gr);

	int a;
	for(int i = 0; i < sx.size(); i++)
	{
		for(int j = 0; j < sy.size(); j++)
		{
			a = check_edge(sx.at(i).e1, sy.at(j).e1); if(a >= 0) add_edge(i * 2 + 0, n + j * 2 + 0, gr);
			a = check_edge(sx.at(i).e1, sy.at(j).e2); if(a >= 0) add_edge(i * 2 + 0, n + j * 2 + 1, gr);
			a = check_edge(sx.at(i).e2, sy.at(j).e1); if(a >= 0) add_edge(i * 2 + 1, n + j * 2 + 0, gr);
			a = check_edge(sx.at(i).e2, sy.at(j).e2); if(a >= 0) add_edge(i * 2 + 1, n + j * 2 + 1, gr);
		}
	}

	a2b.clear();
	b2a.clear();
	for(int i = 0; i < 2 * n; i++) a2b.push_back(i);
	for(int i = 0; i < 2 * n; i++) b2a.push_back(i);

	partners.resize(2 * n);
	for(int i = 0; i < n; i++) partners[i * 2] = i * 2 + 1;
	for(int i = 0; i < n; i++) partners[i * 2 + 1] = i * 2;

	build_gene_list();
	return 0;
}

int distoptimizer::simplify_graph()
{
	int n = num_vertices(gr);
	vector<int> v = partners;

	out_edge_iterator ei1, ei2;
	out_edge_iterator ei3, ei4;
	while(true)
	{
		bool found = false;
		for(int i = 0; i < n; i++)
		{
			if(v[i] == -1) continue;

			tie(ei1, ei2) = out_edges(i, gr);
			tie(ei3, ei4) = out_edges(v[i], gr);

			int d1 = distance(ei1, ei2);
			int d2 = distance(ei3, ei4);

			if(d1 > 1) continue;
			if(d2 > 1) continue;
			found = true;

			int t1 = target(*ei1, gr);
			int t2 = target(*ei3, gr);
			assert(v[t1] >= 0);
			assert(v[t2] >= 0);

			if(v[t1] == t2)
			{
				assert(v[t2] == t1);
				simplified_cycles++;
			}
			else
			{
				assert(v[t2] != t1);	
			}

			add_extremity_pair(b2a[i], b2a[t1]);
			add_extremity_pair(b2a[v[i]], b2a[t2]);

			v[v[t1]] = v[t2];
			v[v[t2]] = v[t1];
			v[i] = v[v[i]] = -1;
			v[t1] = v[t2] = -1;
		}
		if(found == false) break;
	}

	ugraph sg;
	vector<int> c2b;
	vector<int> b2c;

	int cnt = 0;
	for(int i = 0; i < n; i++)
	{
		if(v[i] != -1)
		{
			add_vertex(sg);
			c2b.push_back(i);
			b2c.push_back(cnt);
			cnt++;
		}
		else
		{
			b2c.push_back(-1);
		}
	}
	assert(cnt == num_vertices(sg));
	assert(c2b.size() == num_vertices(sg));
	assert(b2c.size() == num_vertices(gr));

	for(int i = 0; i < n; i++)
	{
		if(b2c[i] == -1) continue;
		for(tie(ei1, ei2) = out_edges(i, gr); ei1 != ei2; ei1++)
		{
			int s = source(*ei1, gr);
			int t = target(*ei1, gr);
			assert(s == i);
			//if(s > t) continue;
			assert(b2c[s] >= 0);
			assert(b2c[t] >= 0);
			add_edge(b2c[s], b2c[t], sg);
		}
	}

	partners.clear();
	partners.assign(num_vertices(sg), -1);
	for(int i = 0; i < n; i++)
	{
		if(v[i] == -1) continue;
		partners[b2c[i]] = b2c[v[i]];
	}

	vector<PII> list;
	for(int i = 0; i < gene_list.size(); i++)
	{
		int x1 = gene_list[i].first.first;
		int x2 = gene_list[i].first.second;
		int g = gene_list[i].second;
		
		if(b2c[x1] == -1) continue;
		if(b2c[x2] == -1) continue;
		
		list.push_back(PII(PI(b2c[x1], b2c[x2]), g));
	}

	vector<int> a2c;
	vector<int> c2a;
	for(int i = 0; i < a2b.size(); i++) 
	{
		int x = a2b[i];
		if(x == -1) a2c.push_back(-1);
		else a2c.push_back(b2c[x]);
	}
	for(int i = 0; i < c2b.size(); i++)
	{
		int x = c2b[i];
		assert(x != -1);
		c2a.push_back(b2a[x]);
	}

	a2b = a2c;
	b2a = c2a;
	gr = sg;
	gene_list = list;

	build_gene_list();

	return 0;
}

int distoptimizer::build_gene_list()
{
	map<gene*, PI> m;
	map<gene*, PI>::iterator it;
	for(int i = 0; i < sx.size(); i++)
	{
		it = m.find(sx.at(i).e1.g);
		if(it == m.end()) 
		{
			if(sx.at(i).e1.b == true) m.insert(pair<gene*, PI>(sx.at(i).e1.g, PI(i * 2, -1)));
			else m.insert(pair<gene*, PI>(sx.at(i).e1.g, PI(-1, i * 2)));
		}
		else
		{
			if(it->second.second == -1) it->second.second = i * 2;
			else if(it->second.first == -1) it->second.first = i * 2;
			else assert(false);
		}

		it = m.find(sx.at(i).e2.g);
		if(it == m.end()) 
		{
			if(sx.at(i).e2.b == true) m.insert(pair<gene*, PI>(sx.at(i).e2.g, PI(i * 2 + 1, -1)));
			else m.insert(pair<gene*, PI>(sx.at(i).e2.g, PI(-1, i * 2 + 1)));
		}
		else
		{
			if(it->second.second == -1) it->second.second = i * 2 + 1;
			else if(it->second.first == -1) it->second.first = i * 2 + 1;
			else assert(false);
		}
	}

	int n = sx.size() * 2;
	for(int i = 0; i < sy.size(); i++)
	{
		it = m.find(sy.at(i).e1.g);
		if(it == m.end()) 
		{
			if(sy.at(i).e1.b == true) m.insert(pair<gene*, PI>(sy.at(i).e1.g, PI(n + i * 2, -1)));
			else m.insert(pair<gene*, PI>(sy.at(i).e1.g, PI(-1, n + i * 2)));
		}
		else
		{
			if(it->second.second == -1) it->second.second = n + i * 2;
			else if(it->second.first == -1) it->second.first = n + i * 2;
			else assert(false);
		}

		it = m.find(sy.at(i).e2.g);
		if(it == m.end()) 
		{
			if(sy.at(i).e2.b == true) m.insert(pair<gene*, PI>(sy.at(i).e2.g, PI(n + i * 2 + 1, -1)));
			else m.insert(pair<gene*, PI>(sy.at(i).e2.g, PI(-1, n + i * 2 + 1)));
		}
		else
		{
			if(it->second.second == -1) it->second.second = i * 2 + 1 + n;
			else if(it->second.first == -1) it->second.first = i * 2 + 1 + n;
			else assert(false);
		}
	}

	gene_list.clear();
	for(it = m.begin(); it != m.end(); it++)
	{
		if(it->second.second == -1)
		{
			assert(it->first->x == 0);
			if(a2b[it->second.first] == -1) continue;
			PI pi(a2b[it->second.first], a2b[it->second.first]);
			gene_list.push_back(PII(pi, (int)fabs(it->first->x)));
		}
		else if(it->second.first == -1)
		{
			assert(it->first->x == 0);
			if(a2b[it->second.second] == -1) continue;
			PI pi(a2b[it->second.second], a2b[it->second.second]);
			gene_list.push_back(PII(pi, (int)fabs(it->first->x)));
		}
		else
		{
			if(a2b[it->second.first] == -1) continue;
			if(a2b[it->second.second] == -1) continue;
			PI pi(a2b[it->second.first], a2b[it->second.second]);
			gene_list.push_back(PII(pi, (int)fabs(it->first->x)));
		}
	}
	return 0;
}

int distoptimizer::check_components()
{
	ugraph g(gr);

	out_edge_iterator ei1, ei2, ei3, ei4;
	for(int i = 0; i < gene_list.size(); i++)
	{
		int x = gene_list[i].first.first;
		int y = gene_list[i].first.second;
		if(x == y) continue;

		tie(ei1, ei2) = out_edges(x, g);
		tie(ei3, ei4) = out_edges(y, g);

		int d1 = distance(ei1, ei2);
		int d2 = distance(ei3, ei4);

		assert(d1 == d2);
		if(d1 == 1) continue;

		//if(d1 != 1) add_edge(x, y, g);
	}

	for(int i = 0; i < partners.size(); i++)
	{
		if(i >= partners[i]) continue;
		add_edge(i, partners[i], g);
	}

	int n = num_vertices(g);
	vector<int> v(n);
	int c = connected_components(g, &v[0]);

	printf("total %d connected components\n", c);

	vector<int> vv;
	vv.assign(c, 0);
	for(int i = 0; i < v.size(); i++) vv[v[i]]++;
	
	for(int i = 0; i < c; i++)
	{
		if(vv[i] == 4) continue;
		printf("component %4d has %4d vertices.\n", i + 1, vv[i]);
	}

	return 0;
}

int distoptimizer::print_breakpoint_graph()
{
	for(int i = 0; i < sx.size() * 2; i++)
	{
		int a = 0, b = 0;
		out_edge_iterator it1, it2;
		for(tie(it1, it2) = out_edges(i, bg); it1 != it2; it1++)
		{
			int x = target(*it1, bg);
			if(x < sx.size() * 2) a++;
			else
			{
				b++;
				printf("edge %2d <-> %2d\n", i, x - sx.size() * 2);
			}
		}
		assert(a == 1 && b == 1);
	}
	return 0;
}
