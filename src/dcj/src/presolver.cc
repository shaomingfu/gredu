#include "presolver.h"
#include "mygraph.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <boost/graph/connected_components.hpp>


using namespace std;

presolver::presolver(ugraph & _gr, const vector<PII> & _gene_list, const vector<int> & _partners)
	: gr(_gr), gene_list(_gene_list), partners(_partners)
{
	//printf("build neighbor indices ...\n");
	build_neighbor_indices();

	//printf("build gene indices ...\n");
	build_gene_indices();

	//printf("build complements ...\n");
	build_complements();

	build_identical_graph();
}

presolver::~presolver()
{
}

bool presolver::simplify(int level)
{
	bool result = false;
	while(true)
	{
		bool b1 = simplify_graph_by_half_fixed_length2();
		//bool b2 = simplify_graph_by_supported_length2(2);
		//bool b3 = simplify_graph_by_checked_length2(2);
		bool b2 = false;
		bool b3 = false;
		bool b4 = false;
		bool b5 = false;
		bool b6 = false;
		//bool b1 = simplify_graph_by_length2();
		//bool b2 = simplify_graph_by_length4();
		//bool b3 = simplify_graph_by_length6();
		//bool b4 = simplify_graph_by_components();
		//bool b3 = simplify_graph_by_pure_length4();

		if(b1 || b2 || b3 || b4 || b5 || b6) result = true;

		//printf("simplify = (%d,%d,%d,%d,%d,%d)\n", b1,b2,b3,b4,b5,b6);
		if(b1 == false && b2 == false && b3 == false && b4 == false && b5 == false && b6 == false) break;
	}

	return result;
}

int presolver::build_neighbor_indices()
{
	int n = num_vertices(gr) / 2;
	neighbors.clear();
	nb_indices.clear();
	nb_indices.assign(2 * n, -1);
	for(int i = 0; i < 2 * n; i++)
	{
		if(nb_indices[i] != -1) continue;
		if(i < partners[i]) neighbors.push_back(PI(i, partners[i]));
		else neighbors.push_back(PI(partners[i], i));
		nb_indices[i] = neighbors.size() - 1;
		nb_indices[partners[i]] = neighbors.size() - 1;
	}
	assert(neighbors.size() == n);
	return 0;
}

int presolver::build_gene_indices()
{
	gn_indices.assign(num_vertices(gr), -1);
	for(int i = 0; i < gene_list.size(); i++)
	{
		gn_indices[gene_list[i].first.first] = i;
		gn_indices[gene_list[i].first.second] = i;

		//printf("%4d -> %4d, gene = %4d\n", gene_list[i].first.first, i, gene_list[i].second);
		//printf("%4d -> %4d, gene = %4d\n", gene_list[i].first.second, i, gene_list[i].second);
	}
	return 0;
}

int presolver::build_complements()
{
	complements.assign(num_vertices(gr), -1);
	for(int i = 0; i < gene_list.size(); i++)
	{
		complements[gene_list[i].first.first] = gene_list[i].first.second;
		complements[gene_list[i].first.second] = gene_list[i].first.first;
	}
	return 0;
}


int presolver::build_identical_graph()
{
	ig.clear();

	assert(neighbors.size() % 2 == 0);
	int n = neighbors.size() / 2;

	for(int i = 0; i < 2 * n; i++) add_vertex(ig);

	for(int i = 0; i < n; i++)
	{
		for(int j = n; j < 2 * n; j++)
		{
			int b = check_identity(i, j);
			if(b <= 0) continue;
			add_edge(i, j, ig);
		}
	}
	return 0;
}

int presolver::analyze_components(const ugraph & g)
{
	vector<int> v( num_vertices(g) );
	int c = connected_components(g, &v[0]);

	vector< vector<int> > vv;
	vv.resize(c);
	for(int i = 0; i < v.size(); i++)
	{
		vv[v[i]].push_back(i);
	}

	int max_size = 0;
	for(int i = 0; i < c; i++)
	{
		if(vv[i].size() > max_size) max_size = vv[i].size();
	}

	vector<int> cc;
	cc.assign(max_size, 0);
	for(int i = 0; i < c; i++)
	{
		cc[vv[i].size() - 1]++;
	}

	for(int i = 0; i < max_size; i++)
	{
		printf("%4d components have %4d adjacencies\n", cc[i], i + 1);
	}

	return 0;
}
bool presolver::simplify_graph_by_fixing()
{
	int sum = 0;
	out_edge_iterator ei1, ei2;
	for(int i = 0; i < num_vertices(gr); i++)
	{
		if(out_degree(i, gr) != 1) continue;
		tie(ei1, ei2) = out_edges(i, gr);
		assert(distance(ei1, ei2) == 1);
		int s = source(*ei1, gr);
		int t = target(*ei1, gr);
		assert(s == i);
		//assert(out_degree(t, gr) == 1);

		if(is_fixed(s, t)) continue;

		sum++;
		if(s < t) fix_extremity_pair(s, t);
		if(t < s) fix_extremity_pair(t, s);
	}
	if(sum >= 1) 
	{
		//printf("fixing %5d pairs\n", sum);
		return true;
	}
	else return false;
}

bool presolver::simplify_graph_by_components()
{
	// TODO, this simplification is dangerous
	bool result = false;
	while(true)
	{
		vector<int> v(neighbors.size());
		int c = connected_components(ig, &v[0]);

		vector< vector<int> > vv;
		vv.resize(c);
		for(int i = 0; i < v.size(); i++) vv[v[i]].push_back(i);

		bool found = false;
		for(int i = 0; i < c; i++)
		{
			bool b = fix_component(vv[i]);
			if(b == true)
			{
				found = true;
				result = true;
				break;
			}
		}

		if(found == false) break;
	}
	return result;
}

bool presolver::remove_extra_edge(int x, const vector<int> & vy)
{
	out_edge_iterator ei1, ei2;
	assert(vy.size() > 0);
	for(tie(ei1, ei2) = out_edges(x, gr); ei1 != ei2; ei1++)
	{
		int s = source(*ei1, gr);
		int t = target(*ei1, gr);
		assert(s == x);

		bool b = false;
		for(int k = 0; k < vy.size(); k++)
		{
			int y1 = neighbors[vy[k]].first;
			int y2 = neighbors[vy[k]].second;
			//assert(edge(x, y1, gr).second || edge(x, y2, gr).second);

			if(t == y1 || t == y2)
			{
				b = true;
				break;
			}
		}
		if(b == true) continue;
		if(out_degree(t, gr) == 1)
		{
			fix_extremity_pair(x, t);
			return true;
		}
		else
		{
			remove_extremity_pair(x, t);
			return true;
		}
	}
	return false;
}

bool presolver::fix_component(const vector<int> & v)
{
	int n = neighbors.size() / 2;
	vector<int> vx;
	vector<int> vy;
	for(int i = 0; i < v.size(); i++)
	{
		if(v[i] <= n) vx.push_back(v[i]);
		else vy.push_back(v[i]);
	}

	if(vx.size() == 1 && vy.size() >= 2)
	{
		int x1 = neighbors[vx[0]].first;
		int x2 = neighbors[vx[0]].second;
		if(remove_extra_edge(x1, vy) == true) return true;
		if(remove_extra_edge(x2, vy) == true) return true;
	}
	else if(vx.size() >= 2 && vy.size() == 1)
	{
		int y1 = neighbors[vy[0]].first;
		int y2 = neighbors[vy[0]].second;
		if(remove_extra_edge(y1, vx) == true) return true;
		if(remove_extra_edge(y2, vx) == true) return true;
	}
	else if(vx.size() == 2 && vy.size() == 2)
	{
		if(edge(vx[0], vy[0], ig).second == false)
		{
			assert(edge(vx[0], vy[1], ig).second);
			assert(edge(vx[1], vy[0], ig).second);
			assert(edge(vx[1], vy[1], ig).second);
			fix_adjacency_pair(vx[0], vy[1]);
			//fix_adjacency_pair(vx[1], vy[0]);
			return true;
		}
		else if(edge(vx[0], vy[1], ig).second == false)
		{
			assert(edge(vx[0], vy[0], ig).second);
			assert(edge(vx[1], vy[1], ig).second);
			assert(edge(vx[1], vy[0], ig).second);
			fix_adjacency_pair(vx[0], vy[0]);
			//fix_adjacency_pair(vx[1], vy[1]);
			return true;
		}
		else if(edge(vx[1], vy[1], ig).second == false)
		{
			assert(edge(vx[0], vy[0], ig).second);
			assert(edge(vx[1], vy[0], ig).second);
			assert(edge(vx[0], vy[1], ig).second);
			fix_adjacency_pair(vx[0], vy[1]);
			//fix_adjacency_pair(vx[1], vy[0]);
			return true;
		}
		else if(edge(vx[1], vy[0], ig).second == false)
		{
			assert(edge(vx[0], vy[0], ig).second);
			assert(edge(vx[1], vy[1], ig).second);
			assert(edge(vx[0], vy[1], ig).second);
			fix_adjacency_pair(vx[0], vy[0]);
			//fix_adjacency_pair(vx[1], vy[1]);
			return true;
		}

		for(int i = 0; i < 2; i++)
		{
			int x1 = neighbors[vx[i]].first;
			int x2 = neighbors[vx[i]].second;
			if(remove_extra_edge(x1, vy) == true) return true;
			if(remove_extra_edge(x2, vy) == true) return true;
		}

		for(int i = 0; i < 2; i++)
		{
			int y1 = neighbors[vy[i]].first;
			int y2 = neighbors[vy[i]].second;
			if(remove_extra_edge(y1, vx) == true) return true;
			if(remove_extra_edge(y2, vx) == true) return true;
		}
	}

	return false;
}

int presolver::remove_extremity_pair(int x, int y)
{
	assert(edge(x, y, gr).second);
	assert(out_degree(x, gr) >= 2);
	assert(out_degree(y, gr) >= 2);

	int cx = complements[x];
	int cy = complements[y];
	assert(cx != -1 && cy != -1);
	assert(out_degree(cx, gr) == out_degree(x, gr));
	assert(out_degree(cy, gr) == out_degree(y, gr));

	assert(edge(x, y, gr).second);
	assert(edge(cx, cy, gr).second);

	remove_edge(x, y, gr);
	if(x != cx) remove_edge(cx, cy, gr);

	int nx = nb_indices[x];
	int ny = nb_indices[y];
	if(edge(nx, ny, ig).second) remove_edge(nx, ny, ig);

	nx = nb_indices[cx];
	ny = nb_indices[cy];
	if(edge(nx, ny, ig).second) remove_edge(nx, ny, ig);
	return 0;
}

bool presolver::simplify_graph_by_half_fixed_length2()
{
	bool result = false;
	while(true)
	{
		bool found = false;
		edge_iterator ei1, ei2;
		int x = -1, y = -1;
		for(tie(ei1, ei2) = edges(ig); ei1 != ei2; ei1++)
		{
			x = source(*ei1, ig);
			y = target(*ei1, ig);
			if(out_degree(x, ig) > 1) continue;
			if(out_degree(y, ig) > 1) continue;

			int x1 = neighbors[x].first;
			int x2 = neighbors[x].second;
			int y1 = neighbors[y].first;
			int y2 = neighbors[y].second;

			if(out_degree(x1, gr) > 1 && out_degree(x2, gr) > 1) continue;
			if(out_degree(y1, gr) > 1 && out_degree(y2, gr) > 1) continue;

			found = true;
			result = true;
			break;
		}
		if(found == false) break;
		assert(x != -1 && y != -1);

		fix_adjacency_pair(x, y);

		clear_vertex(x, ig);
		clear_vertex(y, ig);
	}
	return result;
}

bool presolver::simplify_graph_by_checked_length2(int level)
{
	bool result = false;
	while(true)
	{
		bool found = false;
		edge_iterator ei1, ei2;
		int x = -1, y = -1;
		for(tie(ei1, ei2) = edges(ig); ei1 != ei2; ei1++)
		{
			x = source(*ei1, ig);
			y = target(*ei1, ig);
			//if(out_degree(x, ig) > 1) continue;
			//if(out_degree(y, ig) > 1) continue;

			int x1 = neighbors[x].first;
			int x2 = neighbors[x].second;
			int y1, y2;
			int b = check_identity(x, y);
			if(b == 1)
			{
				y1 = neighbors[y].first;
				y2 = neighbors[y].second;
			}
			else if(b == 2)
			{
				y2 = neighbors[y].first;
				y1 = neighbors[y].second;
			}
			else assert(false);

			assert(out_degree(x1, gr) == out_degree(y1, gr));
			assert(out_degree(x2, gr) == out_degree(y2, gr));

			PI m1(0, 0), m2(0, 0);
			if(out_degree(x1, gr) <= 6) m1 = calc_max_identity_degree(x1, y1);
			if(out_degree(x2, gr) <= 6) m2 = calc_max_identity_degree(x2, y2);

			//printf("(%5d, %5d), degree = (%d, %d), max identity degree = (%d, %d)\n", x1, y1, out_degree(x1, gr), out_degree(y1, gr), m1.first, m1.second);
			//printf("(%5d, %5d), degree = (%d, %d), max identity degree = (%d, %d)\n\n", x2, y2, out_degree(x2, gr), out_degree(y2, gr), m2.first, m2.second);

			int d1 = m1.first - m1.second;
			int d2 = m2.first - m2.second;

			if(level == 2 && (d1 >= 2 && d2 >= 0 || d1 >= 0 && d2 >= 2))
			{
				found = true;
				result = true;
				break;
			}

			if(level == 1 && d1 >= 1 && d2 >= 1)
			{
				found = true;
				result = true;
				break;
			}
		}

		if(found == false) break;
		assert(x != -1 && y != -1);

		fix_adjacency_pair(x, y);

		clear_vertex(x, ig);
		clear_vertex(y, ig);

		simplify_graph_by_half_fixed_length2();
	}
	return result;
}

PI presolver::calc_identity_degree(int x1, int y1)
{
	out_edge_iterator ei1, ei2;
	tie(ei1, ei2) = out_edges(x1, gr);
	assert(distance(ei1, ei2) == 2);
	int y2 = target(*ei1, gr);
	if(y2 == y1) y2 = target(*(++ei1), gr);
	else assert(target(*(++ei1), gr) == y1);

	tie(ei1, ei2) = out_edges(y1, gr);
	assert(distance(ei1, ei2) == 2);
	int x2 = target(*ei1, gr);
	if(x2 == x1) x2 = target(*(++ei1), gr);
	else assert(target(*(++ei1), gr) == x1);

	int c1 = complements[x1];
	int d1 = complements[y1];
	assert(c1 != -1 && d1 != -1);
	assert(edge(c1, d1, gr).second == true);

	int c2 = complements[x2];
	int d2 = complements[y2];
	assert(c2 != -1 && d2 != -1);
	assert(edge(c2, d2, gr).second == true);

	assert(out_degree(x1, gr) == 2);
	assert(out_degree(x2, gr) == 2);
	assert(out_degree(y1, gr) == 2);
	assert(out_degree(y2, gr) == 2);
	assert(out_degree(c1, gr) == 2);
	assert(out_degree(c2, gr) == 2);
	assert(out_degree(d1, gr) == 2);
	assert(out_degree(d2, gr) == 2);

	assert(complements[x1] != -1);
	assert(complements[x2] != -1);
	assert(complements[y1] != -1);
	assert(complements[y2] != -1);
	assert(complements[c1] != -1);
	assert(complements[c2] != -1);
	assert(complements[d1] != -1);
	assert(complements[d2] != -1);

	int p = 0, q = 0;
	if(edge(nb_indices[x1], nb_indices[y1], ig).second) p++;
	if(edge(nb_indices[x2], nb_indices[y2], ig).second) p++;
	if(edge(nb_indices[c1], nb_indices[d1], ig).second) p++;
	if(edge(nb_indices[c2], nb_indices[d2], ig).second) p++;

	if(edge(nb_indices[x1], nb_indices[y2], ig).second) q++;
	if(edge(nb_indices[x2], nb_indices[y1], ig).second) q++;
	if(edge(nb_indices[c1], nb_indices[d2], ig).second) q++;
	if(edge(nb_indices[c2], nb_indices[d1], ig).second) q++;

	/*
	if(edge(complements[x1], complements[y1], gr).second == true) p++;
	if(edge(complements[x2], complements[y2], gr).second == true) p++;
	if(edge(complements[c1], complements[d1], gr).second == true) p++;
	if(edge(complements[c2], complements[d2], gr).second == true) p++;
	if(edge(complements[x1], complements[y2], gr).second == true) q++;
	if(edge(complements[x2], complements[y1], gr).second == true) q++;
	if(edge(complements[c1], complements[d2], gr).second == true) q++;
	if(edge(complements[c2], complements[d1], gr).second == true) q++;
	*/

	return PI(p, q);
}

int presolver::build_subgraph(int x, int y, ugraph & g, ugraph & g2, vector<PII> & gl, vector<int> & pt)
{
	map<int, int> m1;	// map from gr indices to g indices
	map<int, int> m2;	// map from g indices to gr indices

	g.clear();
	g2.clear();
	int x1 = neighbors[x].first;
	int x2 = neighbors[x].second;
	int y1 = neighbors[y].first;
	int y2 = neighbors[y].second;

	assert(complements[x1] != -1);
	assert(complements[x2] != -1);
	assert(complements[y1] != -1);
	assert(complements[y2] != -1);

	out_edge_iterator ei1, ei2;
	set<int> sx;
	set<int> sy;
	for(tie(ei1, ei2) = out_edges(y1, gr); ei1 != ei2; ei1++)
	{
		int s = source(*ei1, gr);
		int t = target(*ei1, gr);
		if(sx.find(t) == sx.end()) sx.insert(t);
		int c = complements[t];
		assert(c != -1);
		if(sx.find(c) == sx.end()) sx.insert(c);
	}
	for(tie(ei1, ei2) = out_edges(y2, gr); ei1 != ei2; ei1++)
	{
		int s = source(*ei1, gr);
		int t = target(*ei1, gr);
		if(sx.find(t) == sx.end()) sx.insert(t);
		int c = complements[t];
		assert(c != -1);
		if(sx.find(c) == sx.end()) sx.insert(c);
	}

	for(tie(ei1, ei2) = out_edges(x1, gr); ei1 != ei2; ei1++)
	{
		int s = source(*ei1, gr);
		int t = target(*ei1, gr);
		if(sy.find(t) == sy.end()) sy.insert(t);

		int c = complements[t];
		assert(c != -1);
		if(sy.find(c) == sy.end()) sy.insert(c);
	}
	for(tie(ei1, ei2) = out_edges(x2, gr); ei1 != ei2; ei1++)
	{
		int s = source(*ei1, gr);
		int t = target(*ei1, gr);
		if(sy.find(t) == sy.end()) sy.insert(t);

		int c = complements[t];
		assert(c != -1);
		if(sy.find(c) == sy.end()) sy.insert(c);
	}

	set<int>::iterator it;
	for(it = sx.begin(); it != sx.end(); it++)
	{
		int a = *it;
		if(m1.find(a) == m1.end())
		{
			m1.insert(PI(a, num_vertices(g)));
			m2.insert(PI(num_vertices(g), a));
			add_vertex(g);
		}
		int p = partners[a];
		if(m1.find(p) == m1.end())
		{
			m1.insert(PI(p, num_vertices(g)));
			m2.insert(PI(num_vertices(g), p));
			add_vertex(g);
		}
	}

	for(it = sy.begin(); it != sy.end(); it++)
	{
		int a = *it;
		if(m1.find(a) == m1.end())
		{
			m1.insert(PI(a, num_vertices(g)));
			m2.insert(PI(num_vertices(g), a));
			add_vertex(g);
		}
		int p = partners[a];
		if(m1.find(p) == m1.end())
		{
			m1.insert(PI(p, num_vertices(g)));
			m2.insert(PI(num_vertices(g), p));
			add_vertex(g);
		}
	}

	assert(m1.find(x1) != m1.end());
	assert(m1.find(x2) != m1.end());
	assert(m1.find(y1) != m1.end());
	assert(m1.find(y2) != m1.end());

	int n = num_vertices(g);
	for(int i = 0; i < n; i++)
	{
		for(int j = i + 1; j < n; j++)
		{
			int a = m2[i];
			int b = m2[j];
			if(edge(a, b, gr).second == false) continue;
			add_edge(i, j, g);
		}
	}

	g2 = g;
	int b = check_identity(x, y);
	if(b == 1)
	{
		int sx1 = m1[x1];
		int sy1 = m1[y1];
		int tx1 = m1[complements[x1]];
		int ty1 = m1[complements[y1]];
		if(edge(sx1, sy1, g2).second) remove_edge(sx1, sy1, g2);
		if(edge(tx1, ty1, g2).second) remove_edge(tx1, ty1, g2);
		int sx2 = m1[x2];
		int sy2 = m1[y2];
		int tx2 = m1[complements[x2]];
		int ty2 = m1[complements[y2]];
		if(edge(sx2, sy2, g2).second) remove_edge(sx2, sy2, g2);
		if(edge(tx2, ty2, g2).second) remove_edge(tx2, ty2, g2);
	}
	else if(b == 2)
	{
		int sx1 = m1[x1];
		int sy1 = m1[y2];
		int tx1 = m1[complements[x1]];
		int ty1 = m1[complements[y2]];
		if(edge(sx1, sy1, g2).second) remove_edge(sx1, sy1, g2);
		if(edge(tx1, ty1, g2).second) remove_edge(tx1, ty1, g2);
		int sx2 = m1[x2];
		int sy2 = m1[y1];
		int tx2 = m1[complements[x2]];
		int ty2 = m1[complements[y1]];
		if(edge(sx2, sy2, g2).second) remove_edge(sx2, sy2, g2);
		if(edge(tx2, ty2, g2).second) remove_edge(tx2, ty2, g2);
	}
	else assert(false);

	assert(n % 2 == 0);
	pt.resize(n);
	for(int i = 0; i < n; i++)
	{
		assert(m2.find(i) != m2.end());
		int a = m2[i];
		int b = partners[a];
		assert(m1.find(b) != m1.end());
		int c = m1[b];
		assert(c < n);
		pt[i] = c;
	}

	for(int i = 0; i < n; i++)
	{
		int p = pt[i];
		assert(p < n);
		assert(pt[p] == i);
	}

	set<int> gs;
	for(int i = 0; i < n; i++)
	{
		int a = m2[i];
		int gn = gn_indices[a];
		if(gn < 0) continue;

		if(gs.find(gn) != gs.end()) continue;
		gs.insert(gn);

		PII pii = gene_list[gn];
		if(m1.find(pii.first.first) == m1.end()) continue;
		if(m1.find(pii.first.second) == m1.end()) continue;
		pii.first.first = m1[pii.first.first];
		pii.first.second = m1[pii.first.second];
		gl.push_back(pii);
	}

	// print
	
	
	printf("total %5d vertices, %5d edges, %5d pairs\n", (int)num_vertices(g), (int)num_edges(g), (int)gl.size());
	printf("total %5d vertices, %5d edges, %5d pairs\n", (int)num_vertices(g2), (int)num_edges(g2), (int)gl.size());
	return 0;
}

PI presolver::calc_max_identity_degree(int x, int y)
{
	assert(out_degree(x, gr) == out_degree(y, gr));

	set<int> xv;
	set<int> yv;
	out_edge_iterator ei1, ei2;
	assert(edge(x, y, gr).second);
	for(tie(ei1, ei2) = out_edges(x, gr); ei1 != ei2; ei1++)
	{
		int s = source(*ei1, gr);
		int t = target(*ei1, gr);
		assert(s == x);
		if(t == y) continue;
		yv.insert(t);
	}
	for(tie(ei1, ei2) = out_edges(y, gr); ei1 != ei2; ei1++)
	{
		int s = source(*ei1, gr);
		int t = target(*ei1, gr);
		assert(s == y);
		if(t == x) continue;
		xv.insert(t);
	}
	int p = calc_max_identity_degree(x, y, xv, yv);

	set<int>::iterator it;
	int q = 0;
	for(it = yv.begin(); it != yv.end(); it++)
	{
		int yy = *it;
		set<int> yyv = yv;
		yyv.insert(y);
		yyv.erase(yy);
		assert(yy != y);
		int a = calc_max_identity_degree(x, yy, xv, yyv);
		if(a > q) q = a;
	}

	return PI(p, q);
}

int presolver::calc_max_identity_degree(int x, int y, const set<int> & xv, const set<int> & yv)
{
	assert(xv.size() == yv.size());
	assert(xv.find(x) == xv.end());
	assert(yv.find(y) == yv.end());

	set<int>::const_iterator it;
	for(it = xv.begin(); it != xv.end(); it++) assert(edge(*it, y, gr).second);
	for(it = yv.begin(); it != yv.end(); it++) assert(edge(*it, x, gr).second);

	int p = 0;
	assert(edge(x, y, gr).second);
	
	if( edge(nb_indices[x], nb_indices[y], ig).second ) p++;

	int cx = complements[x];
	int cy = complements[y];
	if(cx == -1) assert(cy == -1);
	else assert(cy != -1);

	if(cx >= 0 && edge(nb_indices[cx], nb_indices[cy], ig).second) p++;

	int max = 0;
	int xx = *(xv.begin());
	set<int> xxv = xv;
	xxv.erase(xx);
	for(it = yv.begin(); it != yv.end(); it++)
	{
		int yy = *it;
		set<int> yyv = yv;
		yyv.erase(yy);

		int a = calc_max_identity_degree(xx, yy, xxv, yyv);
		if(a > max) max = a;
	}
	return max + p;
}

bool presolver::simplify_graph_by_length2()
{
	bool result = false;
	while(true)
	{
		bool found = false;
		edge_iterator ei1, ei2;
		int x = -1, y = -1;
		for(tie(ei1, ei2) = edges(ig); ei1 != ei2; ei1++)
		{
			x = source(*ei1, ig);
			y = target(*ei1, ig);
			if(out_degree(x, ig) > 1) continue;
			if(out_degree(y, ig) > 1) continue;
			found = true;
			result = true;
			break;
		}
		if(found == false) break;
		assert(x != -1 && y != -1);

		fix_adjacency_pair(x, y);

		clear_vertex(x, ig);
		clear_vertex(y, ig);
	}
	return result;
}

bool presolver::simplify_graph_by_length4()
{
	bool result = false;
	while(true)
	{
		bool found = false;
		for(int i = 0; i < num_vertices(gr) / 2; i++) 
		{
			int d = out_degree(i, gr);
			assert(d >= 1);
			if(d > 1) continue;

			bool b = fix_length4_cycle(i);
			if(b == true)
			{
				found = true;
				result = true;
				break;
			}
		}
		if(found == false) break;
	}
	return result;
}

bool presolver::simplify_graph_by_pure_length4()
{
	bool result = false;
	while(true)
	{
		bool found = false;
		for(int i = 0; i < num_vertices(gr) / 2; i++) 
		{
			bool b = fix_pure_length4_cycle(i);
			if(b == true)
			{
				found = true;
				result = true;
				break;
			}
		}
		if(found == false) break;
	}
	return result;
}

bool presolver::simplify_graph_by_length6()
{
	bool result = false;
	while(true)
	{
		bool found = false;
		for(int i = 0; i < num_vertices(gr) / 2; i++) 
		{
			int d = out_degree(i, gr);
			assert(d >= 1);
			if(d > 1) continue;

			bool b = fix_length6_cycle(i);
			if(b == true)
			{
				found = true;
				result = true;
				break;
			}
		}
		if(found == false) break;
	}
	return result;
}

bool presolver::unique_length4_cycle(int x1)
{
	if(out_degree(nb_indices[x1], ig) >= 1) return false;

	int cnt = 0;
	vector<PI> ps;

	int x2 = partners[x1];
	out_edge_iterator ei1, ei2;
	out_edge_iterator ei3, ei4;
	out_edge_iterator ei5, ei6;
	for(tie(ei1, ei2) = out_edges(x1, gr); ei1 != ei2; ei1++)
	{
		assert(source(*ei1, gr) == x1);
		int y1 = target(*ei1, gr);
		int y3 = partners[y1];
		if(out_degree(nb_indices[y1], ig) >= 1) return false;

		for(tie(ei3, ei4) = out_edges(y3, gr); ei3 != ei4; ei3++)
		{
			assert(source(*ei3, gr) == y3);
			int x3 = target(*ei3, gr);
			int x4 = partners[x3];
			if(out_degree(nb_indices[x3], ig) >= 1) return false;

			if(x3 == x1 || x3 == x2) continue;
			if(x4 == x1 || x4 == x2) continue;
			for(tie(ei5, ei6) = out_edges(x4, gr); ei5 != ei6; ei5++)
			{
				assert(source(*ei5, gr) == x4);
				int y4 = target(*ei5, gr);
				int y2 = partners[y4];
				if(out_degree(nb_indices[y2], ig) >= 1) return false;
				
				if(y2 == y1 || y2 == y3) continue;
				if(y4 == y1 || y4 == y3) continue;
				if(edge(x2, y2, gr).second == true)
				{
					cnt++;
					if(cnt >= 2) return false;
				}
			}
		}
	}
	if(cnt == 1) return true;
	else return false;
}

bool presolver::fix_pure_length4_cycle(int x1)
{
	if(out_degree(nb_indices[x1], ig) >= 1) return false;

	int cnt = 0;
	vector<PI> ps;

	int x2 = partners[x1];
	out_edge_iterator ei1, ei2;
	out_edge_iterator ei3, ei4;
	out_edge_iterator ei5, ei6;
	for(tie(ei1, ei2) = out_edges(x1, gr); ei1 != ei2; ei1++)
	{
		assert(source(*ei1, gr) == x1);
		int y1 = target(*ei1, gr);
		int y3 = partners[y1];
		if(out_degree(nb_indices[y1], ig) >= 1) return false;

		for(tie(ei3, ei4) = out_edges(y3, gr); ei3 != ei4; ei3++)
		{
			assert(source(*ei3, gr) == y3);
			int x3 = target(*ei3, gr);
			int x4 = partners[x3];
			if(out_degree(nb_indices[x3], ig) >= 1) return false;

			if(x3 == x1 || x3 == x2) continue;
			if(x4 == x1 || x4 == x2) continue;
			for(tie(ei5, ei6) = out_edges(x4, gr); ei5 != ei6; ei5++)
			{
				assert(source(*ei5, gr) == x4);
				int y4 = target(*ei5, gr);
				int y2 = partners[y4];
				if(out_degree(nb_indices[y2], ig) >= 1) return false;
				
				if(y2 == y1 || y2 == y3) continue;
				if(y4 == y1 || y4 == y3) continue;
				if(edge(x2, y2, gr).second == true)
				{
					cnt++;
					if(cnt >= 2) return false;

					ps.clear();
					ps.push_back(PI(x1, y1));
					ps.push_back(PI(x2, y2));
					ps.push_back(PI(x3, y3));
					ps.push_back(PI(x4, y4));

					int dx1 = out_degree(x1, gr);
					int dx2 = out_degree(x2, gr);
					int dx3 = out_degree(x3, gr);
					int dx4 = out_degree(x4, gr);
					int dy1 = out_degree(y1, gr);
					int dy2 = out_degree(y2, gr);
					int dy3 = out_degree(y3, gr);
					int dy4 = out_degree(y4, gr);

					//printf("length4 cycle: (%4d,%4d,%4d,%4d,%4d,%4d,%4d,%4d): ", x1,x2,x3,x4,y1,y2,y3,y4);
					//printf("(%4d,%4d,%4d,%4d,%4d,%4d,%4d,%4d)\n", dx1,dx2,dx3,dx4,dy1,dy2,dy3,dy4);
				}
			}
		}
	}
	if(ps.size() == 0) return false;
	assert(ps.size() == 4);

	assert(unique_length4_cycle(ps[0].first));
	if(unique_length4_cycle(ps[0].second) == false) return false;
	if(unique_length4_cycle(ps[1].first) == false) return false;
	if(unique_length4_cycle(ps[1].second) == false) return false;
	if(unique_length4_cycle(ps[2].first) == false) return false;
	if(unique_length4_cycle(ps[2].second) == false) return false;
	if(unique_length4_cycle(ps[3].first) == false) return false;
	if(unique_length4_cycle(ps[3].second) == false) return false;

	if(is_fixed(ps[0].first, ps[0].second) == false) { fix_extremity_pair(ps[0].first, ps[0].second);  return true; }
	if(is_fixed(ps[1].first, ps[1].second) == false) { fix_extremity_pair(ps[1].first, ps[1].second);  return true; }
	if(is_fixed(ps[2].first, ps[2].second) == false) { fix_extremity_pair(ps[2].first, ps[2].second);  return true; }
	if(is_fixed(ps[3].first, ps[3].second) == false) { fix_extremity_pair(ps[3].first, ps[3].second);  return true; }
	return false;
}

bool presolver::fix_length4_cycle(int x)
{
	assert(out_degree(x, gr) == 1);
	out_edge_iterator ei1, ei2;
	tie(ei1, ei2) = out_edges(x, gr);
	assert(distance(ei1, ei2) == 1);
	assert(source(*ei1, gr) == x);
	int y = target(*ei1, gr);

	int cx = partners[x];
	int cy = partners[y];
	if(edge(cx, cy, gr).second) return false;
	if(edge(x,  cy, gr).second) return false;
	if(edge(cx,  y, gr).second) return false;

	vector<int> vx;
	vector<int> vy;
	for(tie(ei1, ei2) = out_edges(cx, gr); ei1 != ei2; ei1++)
	{
		assert(source(*ei1, gr) == cx);
		vy.push_back(target(*ei1, gr));
	}
	for(tie(ei1, ei2) = out_edges(cy, gr); ei1 != ei2; ei1++)
	{
		assert(source(*ei1, gr) == cy);
		vx.push_back(target(*ei1, gr));
	}

	set<int> nx;
	set<int> ny;
	for(int i = 0; i < vx.size(); i++)
	{
		for(int j = 0; j < vy.size(); j++)
		{
			cx = partners[vx[i]];
			cy = partners[vy[j]];
			if(edge(cx, cy, gr).second)
			{
				nx.insert(vx[i]);
				ny.insert(vy[j]);
			}
		}
	}

	if(nx.size() == 0 && ny.size() == 0) return false;
	if(nx.size() >= 2 && ny.size() >= 2) return false;

	cx = partners[x];
	cy = partners[y];

	if(nx.size() == 1)
	{
		int xx = *(nx.begin());
		int cxx = partners[xx];
		assert(edge(xx, cy, gr).second);
		//if(pairs.find(PI(xx, cy)) == pairs.end() && pairs.find(PI(cy, xx)) == pairs.end())
		if(is_fixed(xx, cy) == false)
		{
			fix_extremity_pair(xx, cy);
			return true;
		}
	}
	if(ny.size() == 1)
	{
		int yy = *(ny.begin());
		int cyy = partners[yy];
		assert(edge(yy, cx, gr).second);
		//if(pairs.find(PI(yy, cx)) == pairs.end() && pairs.find(PI(cx, yy)) == pairs.end())
		if(is_fixed(cx, yy) == false)
		{
			fix_extremity_pair(cx, yy);
			return true;
		}
	}
	return false;
}


bool presolver::fix_length6_cycle(int x)
{
	assert(x >= 0 && x < num_vertices(gr) / 2);
	assert(out_degree(x, gr) == 1);
	out_edge_iterator ei1, ei2;
	tie(ei1, ei2) = out_edges(x, gr);
	assert(distance(ei1, ei2) == 1);
	assert(source(*ei1, gr) == x);
	int y = target(*ei1, gr);

	int cx = partners[x];
	int cy = partners[y];
	if(edge(cx, cy, gr).second) return false;
	if(edge(x,  cy, gr).second) return false;
	if(edge(cx,  y, gr).second) return false;

	if(out_degree(cx, gr) != 1) return false;
	if(out_degree(cy, gr) != 1) return false;

	tie(ei1, ei2) = out_edges(cx, gr);
	assert(distance(ei1, ei2) == 1);
	assert(source(*ei1, gr) == cx);
	int yy = target(*ei1, gr);
	assert(yy != y);

	tie(ei1, ei2) = out_edges(cy, gr);
	assert(distance(ei1, ei2) == 1);
	assert(source(*ei1, gr) == cy);
	int xx = target(*ei1, gr);
	assert(xx != x);

	int cxx = partners[xx];
	int cyy = partners[yy];
	if(edge(cxx, y, gr).second) return false;
	if(edge(cxx, cy, gr).second) return false;
	if(edge(cxx, yy, gr).second) return false;
	if(edge(cxx, cyy, gr).second) return false;
	if(edge(cyy, x, gr).second) return false;
	if(edge(cyy, cx, gr).second) return false;
	if(edge(cyy, xx, gr).second) return false;
	if(edge(cyy, cxx, gr).second) return false;

	vector<int> vx;
	vector<int> vy;
	for(tie(ei1, ei2) = out_edges(cxx, gr); ei1 != ei2; ei1++)
	{
		assert(source(*ei1, gr) == cxx);
		vy.push_back(target(*ei1, gr));
	}
	for(tie(ei1, ei2) = out_edges(cyy, gr); ei1 != ei2; ei1++)
	{
		assert(source(*ei1, gr) == cyy);
		vx.push_back(target(*ei1, gr));
	}

	set<int> nx;
	set<int> ny;
	for(int i = 0; i < vx.size(); i++)
	{
		for(int j = 0; j < vy.size(); j++)
		{
			int u = partners[vx[i]];
			int v = partners[vy[j]];
			if(edge(u, v, gr).second)
			{
				nx.insert(vx[i]);
				ny.insert(vy[j]);
			}
		}
	}

	//printf("x = %5d, y = %5d, vx = %5d, vy = %5d, nx = %5d, ny = %5d\n", x, y, vx.size(), vy.size(), nx.size(), ny.size());

	if(nx.size() == 0 && ny.size() == 0) return false;
	if(nx.size() >= 2 && ny.size() >= 2) return false;

	if(nx.size() == 1)
	{
		int xxx = *(nx.begin());
		int cxxx = partners[xxx];
		assert(edge(xxx, cyy, gr).second);
		//if(pairs.find(PI(xxx, cyy)) == pairs.end() && pairs.find(PI(cyy, xxx)) == pairs.end())
		if(is_fixed(xxx, cyy) == false)
		{
			fix_extremity_pair(xxx, cyy);
			return true;
		}
	}
	if(ny.size() == 1)
	{
		int yyy = *(ny.begin());
		int cyyy = partners[yyy];
		assert(edge(yyy, cxx, gr).second);
		//if(pairs.find(PI(yyy, cxx)) == pairs.end() && pairs.find(PI(cxx, yyy)) == pairs.end())
		if(is_fixed(cxx, yyy) == false)
		{
			fix_extremity_pair(cxx, yyy);
			return true;
		}
	}
	return false;
}

int presolver::fix_adjacency_pair(int x, int y)
{
	int x1 = neighbors[x].first;
	int x2 = neighbors[x].second;
	int y1 = neighbors[y].first;
	int y2 = neighbors[y].second;

	int b = check_identity(x, y);
	if(b == 1)
	{
		fix_extremity_pair(x1, y1);
		fix_extremity_pair(x2, y2);
	}
	else if(b == 2)
	{
		fix_extremity_pair(x1, y2);
		fix_extremity_pair(x2, y1);
	}
	else assert(false);

	return 0;
}

int presolver::fix_extremity(int x)
{
	assert(out_degree(x, gr) == 1);
	out_edge_iterator ei1, ei2;
	tie(ei1, ei2) = out_edges(x, gr);
	assert(distance(ei1, ei2) == 1);
	assert(source(*ei1, gr) == x);
	fix_extremity_pair(x, target(*ei1, gr));
	return 0;
}

bool presolver::is_fixed(int x, int y)
{
	if(x > y) return is_fixed(y, x);

	int n = num_vertices(gr) / 2;
	assert(x >= 0 && x < n);
	assert(y >= n && y < 2 * n);

	pair<edge_descriptor, bool> p = edge(x, y, gr);
	assert(p.second == true);

	if(out_degree(x, gr) > 1) return false;
	if(out_degree(y, gr) > 1) return false;

	int cx = complements[x];
	int cy = complements[y];

	if(cx == -1) assert(cy == -1);
	else assert(cy != -1);
	if(cx == -1) return true;

	assert(edge(cx, cy, gr).second);
	assert(cx < cy);

	if(out_degree(cx, gr) > 1) return false;
	if(out_degree(cy, gr) > 1) return false;

	return true;
}

int presolver::fix_extremity_pair(int x, int y)
{
	int n = num_vertices(gr) / 2;
	assert(x >= 0 && x < n);
	assert(y >= n && y < 2 * n);

	pair<edge_descriptor, bool> p = edge(x, y, gr);
	assert(p.second == true);

	clear_vertex(x, gr);
	clear_vertex(y, gr);
	add_edge(x, y, gr);

	// update identity graph
	int nx = nb_indices[x];
	int ny = nb_indices[y];
	p = edge(nx, ny, ig);
	clear_vertex(nx, ig);
	clear_vertex(ny, ig);
	if(p.second) add_edge(nx, ny, ig);

	int cx = complements[x];
	int cy = complements[y];

	if(cx == -1) assert(cy == -1);
	else assert(cy != -1);
	if(cx == -1) return 0;

	assert(edge(cx, cy, gr).second);

	assert(cx < cy);

	clear_vertex(cx, gr);
	clear_vertex(cy, gr);
	add_edge(cx, cy, gr);

	// update identity graph
	nx = nb_indices[cx];
	ny = nb_indices[cy];

	p = edge(nx, ny, ig);
	clear_vertex(nx, ig);
	clear_vertex(ny, ig);

	if(p.second) add_edge(nx, ny, ig);

	return 0;
}

int presolver::check_identity(int x, int y)
{
	int n = neighbors.size() / 2;
	//assert(x >= 0 && x < n);
	//assert(y >= n && y < 2 * n);

	int x1 = neighbors[x].first;
	int x2 = neighbors[x].second;
	int y1 = neighbors[y].first;
	int y2 = neighbors[y].second;

	int cx1 = complements[x1];
	int cy1 = complements[y1];

	if(cx1 == x2 && cy1 != y2) return -2;
	if(cx1 != x2 && cy1 == y2) return -2;

	pair<edge_descriptor, bool> p1, p2, p3, p4;
	p1 = edge(x1, y1, gr);
	p2 = edge(x2, y2, gr);
	p3 = edge(x1, y2, gr);
	p4 = edge(x2, y1, gr);
	if(p1.second && p2.second && p3.second && p4.second) return 0;
	if(p1.second && p2.second) return 1;
	if(p3.second && p4.second) return 2;
	return -1;
}

int presolver::print()
{
	// print graph
	printf("================= graph =================\n");
	for(int i = 0; i < num_vertices(gr); i++)
	{
		int g = gn_indices[i];
		char c;
		if(g == -1) c = '-';
		else if(gene_list[g].first.first == i) c = 'H';
		else if(gene_list[g].first.second == i) c = 'T';
		else assert(false);

		int gg = -1;
		if(g != -1) gg = gene_list[g].second;

		printf("vertex %4d [partner = %4d] [complement = %4d] [%4d:%4d:%c]: ", i, partners[i], complements[i], g, gg, c); 

		out_edge_iterator ei1, ei2;
		for(tie(ei1, ei2) = out_edges(i, gr); ei1 != ei2; ei1++)
		{
			int s = source(*ei1, gr);
			int t = target(*ei1, gr);
			assert(s == i);
			printf("%4d, ", t);
		}
		printf("\n");
	}

	printf("================= gene list =================\n");
	for(int i = 0; i < gene_list.size(); i++)
	{
		printf("gene: %4d, (%4d,%4d): %4d\n", i, gene_list[i].first.first, gene_list[i].first.second, gene_list[i].second);
	}

	return 0;
}

int presolver::build_connected_graph()
{
	cg.clear();
	int n = num_vertices(ig);
	for(int i = 0; i < n; i++) add_vertex(cg);

	out_edge_iterator ei1, ei2;
	for(int i = 0; i < n; i++)
	{
		if(out_degree(i, gr) > 1) continue;
		tie(ei1, ei2) = out_edges(i, gr);
		assert(distance(ei1, ei2) == 1);
		assert(source(*ei1, gr) == i);
		int j = target(*ei1, gr);

		int ni = nb_indices[i];
		int nj = nb_indices[j];
		assert(ni < n / 2);
		assert(nj >= n / 2);

		add_edge(ni, nj, cg);
	}

	analyze_components(cg);

	int r = 0;
	while(true)
	{
		printf("round %d merge\n", r++);
		bool b = merge_connected_graph();
		analyze_components(cg);
		if(b == false) break;
	}
	return 0;
}

set<int> presolver::get_linked_adjacencies(int x1)
{
	set<int> s;
	out_edge_iterator ei1, ei2;
	for(tie(ei1, ei2) = out_edges(x1, gr); ei1 != ei2; ei1++)
	{
		int t = target(*ei1, gr);
		int nt = nb_indices[t];
		if(s.find(nt) == s.end()) s.insert(nt);
	}
	return s;
}

bool presolver::merge_connected_graph()
{
	int n = num_vertices(cg);
	vector<int> v(n);
	int c = connected_components(cg, &v[0]);

	out_edge_iterator ei1, ei2;
	for(int i = 0; i < n; i++)
	{
		int ni = nb_indices[i];
		set<int> si = get_linked_adjacencies(i);
		if(si.size() != 2) continue;
		set<int>::iterator it = si.begin();
		int x = *it;
		int y = *(++it);

		assert(ni < n / 2);
		assert(x >= n / 2);
		assert(y >= n / 2);

		int x1 = neighbors[x].first;
		int x2 = neighbors[x].second;
		set<int> sx1 = get_linked_adjacencies(x1);
		set<int> sx2 = get_linked_adjacencies(x2);
		int y1 = neighbors[y].first;
		int y2 = neighbors[y].second;
		set<int> sy1 = get_linked_adjacencies(y1);
		set<int> sy2 = get_linked_adjacencies(y2);

		if(sx1.size() == 2)
		{
			printf("AAAAAAAAAAAA1\n");
			it = sx1.begin();
			int a = *it;
			int b = *(++it);
			printf("%5d: (%5d, %5d, %5d, %5d) : (%5d, %5d, %5d, %5d)\n", ni, a, b, x, y, v[a], v[b], v[x], v[y]);
			if(a == ni && v[b] == v[y] && edge(ni, x, cg).second == false)
			{
				add_edge(ni, x, cg);
				return true;
			}
			else if(b == ni && v[a] == v[y] && edge(ni, x, cg).second == false)
			{
				add_edge(ni, x, cg);
				return true;
			}
		}
		if(sx2.size() == 2)
		{
			printf("AAAAAAAAAAAA2\n");
			it = sx2.begin();
			int a = *it;
			int b = *(++it);
			if(a == ni && v[b] == v[y] && edge(ni, x, cg).second == false)
			{
				add_edge(ni, x, cg);
				return true;
			}
			else if(b == ni && v[a] == v[y] && edge(ni, x, cg).second == false)
			{
				add_edge(ni, x, cg);
				return true;
			}
		}
		if(sy1.size() == 2)
		{
			printf("AAAAAAAAAAAA3\n");
			it = sy1.begin();
			int a = *it;
			int b = *(++it);
			if(a == ni && v[b] == v[x] && edge(ni, y, cg).second == false)
			{
				add_edge(ni, y, cg);
				return true;
			}
			else if(b == ni && v[a] == v[x] && edge(ni, y, cg).second == false)
			{
				add_edge(ni, y, cg);
				return true;
			}
		}
		if(sy2.size() == 2)
		{
			printf("AAAAAAAAAAAA4\n");
			it = sy2.begin();
			int a = *it;
			int b = *(++it);
			if(a == ni && v[b] == v[x] && edge(ni, y, cg).second == false)
			{
				add_edge(ni, y, cg);
				return true;
			}
			else if(b == ni && v[a] == v[x] && edge(ni, y, cg).second == false)
			{
				add_edge(ni, y, cg);
				return true;
			}
		}
	}
	return false;
}


