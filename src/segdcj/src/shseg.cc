#include "draw.h"
#include "shseg.h"
#include <algorithm>
#include <cassert>
#include <cstdio>
#include <cmath>

shseg::shseg()
{}

shseg::shseg(const vector<int> & x, const vector<int> & y)
	: xv(x), yv(y)
{
	assert(xv.size() == yv.size());
}

shseg::shseg(const vector<PI> & v)
{
	xv.clear();
	yv.clear();
	for(int i = 0; i < v.size(); i++)
	{
		xv.push_back(v[i].first);
		yv.push_back(v[i].second);
	}
}

bool shseg::operator< (const shseg & ss) const
{
	if(eb.size() < ss.eb.size()) return true;
	else return false;
}

int shseg::size() const
{
	return xv.size();
}

int shseg::calc_lower_bound()
{
	ugraph gr;
	int n = num_vertices(ar);
	for(int i = 0; i < n; i++) add_vertex(gr);
	for(int i = 0; i < adj.size(); i++)
	{
		int j = adj[i];
		if(j == -1) continue;
		assert(adj[j] == i);
		if(i > j) continue;
		assert(i < n && j < n);
		add_edge(i, j, gr);
	}

	for(set<PI>::iterator it = eb.begin(); it != eb.end(); it++)
	{
		int x = it->first;
		int y = it->second;
		add_edge(ei[x], ei[y], gr);
	}

	for(int i = 0; i < xv.size(); i++)
	{
		int x = xv[i];
		int y = yv[i];
		add_edge(ei[2 * x + 0], ei[2 * y + 0], gr);
		add_edge(ei[2 * x + 1], ei[2 * y + 1], gr);
	}

	for(int i = 0; i < ge.size(); i++)
	{
		int x = ge[i].first;
		int y = ge[i].second;
		int g = ie[x] / 2;
		bool b = false;
		for(int k = 0; k < xv.size(); k++)
		{
			if(xv[k] == g) b = true;
			if(yv[k] == g) b = true;
			if(b == true) break;
		}
		if(b == true) continue;
		add_edge(x, y, gr);
	}

	vector<int> v(num_vertices(gr));
	int c = connected_components(gr, &v[0]);
	return c;
}

int shseg::statistic() const
{
	printf("segment has %2d genes, %4d intra-extremities, %4d extra extremities, %4d intra-edges %4d extra edges\n",
			(int)xv.size(), (int)sa.size(), (int)sb.size(), (int)ea.size(), (int)eb.size());
	return 0;
}

int shseg::print() const
{
	printf("genes:\n");
	assert(ge.size() == gf.size());
	for(int i = 0; i < ge.size(); i++)
	{
		int x = ge[i].first;
		int y = ge[i].second;
		printf("gene %3d = (%3d,%3d), adj = (%3d,%3d), fam = %3d\n", i, x, y, adj[x], adj[y], gf[i]);
	}

	printf("edges:\n");
	edge_iterator ei1, ei2;
	int cnt = 0;
	for(tie(ei1, ei2) = edges(ar); ei1 != ei2; ei1++)
	{
		int s = source(*ei1, ar);
		int t = target(*ei1, ar);
		printf("edge %3d = (%3d,%3d)\n", cnt++, s, t);
	}

	printf("duplicons:\n");
	for(int i = 0; i < vdup.size(); i++)
	{
		printf("duplicon %3d: ", i);
		for(int k = 0; k < vdup[i].size(); k++)
		{
			printf("%3d ", vdup[i][k]);
		}
		printf("\n");
	}
	return 0;
}

bool shseg::contain(int x, int y) const
{
	for(int i = 0; i < xv.size(); i++)
	{
		if(xv[i] == x && yv[i] == y) return true;
	}
	return false;
}

map<int, int> shseg::sort_vertices(int xsize, bool b)
{
	map<int, int> mi;
	map<int, int> im;
	for(map<int, int>::iterator it = ei.begin(); it != ei.end(); it++)
	{
		if(b == true && it->first >= xsize) continue;
		if(b == false && it->first < xsize) continue;
		int k = mi.size();
		mi.insert(PI(it->second, k));
		im.insert(PI(k, it->second));
	}

	ugraph gr;
	for(int i = 0; i < mi.size(); i++) add_vertex(gr);
	for(int i = 0; i < adj.size(); i++)
	{
		if(mi.find(i) == mi.end()) continue;
		int j = adj[i];
		assert(mi.find(j) != mi.end());
		if(i > j) continue;
		assert(edge(mi[i], mi[j], gr).second == false);
		add_edge(mi[i], mi[j], gr);
	}
	for(int i = 0; i < ge.size(); i++)
	{
		int x = ge[i].first;
		int y = ge[i].second;
		if(mi.find(x) == mi.end()) continue;
		assert(mi.find(y) != mi.end());
		if(edge(mi[x], mi[y], gr).second) continue;
		add_edge(mi[x], mi[y], gr);
	}

	for(int i = 0; i < num_vertices(gr); i++)
	{
		assert(out_degree(i, gr) >= 1);
		assert(out_degree(i, gr) <= 2);
	}

	vector<int> v;
	v.assign(mi.size(), -1);
	int index = 0;
	while(true)
	{
		int start = -1;
		for(int i = 0; i < mi.size(); i++)
		{
			if(out_degree(i, gr) == 1 && v[i] == -1)
			{
				start = i;
				break;
			}
		}
		if(start == -1)
		{
			for(int i = 0; i < v.size(); i++)
			{
				if(v[i] != -1) continue;
				assert(out_degree(i, gr) == 2);
				start = i;
				break;
			}
		}
		if(start == -1) break;

		v[start] = index;
		index++;

		int p = start;
		while(true)
		{
			out_edge_iterator ei1, ei2;
			int q = -1;
			for(tie(ei1, ei2) = out_edges(p, gr); ei1 != ei2; ei1++)
			{
				int x = target(*ei1, gr);
				if(v[x] != -1) continue;
				q = x;
				break;
			}
			if(q == -1) break;
			v[q] = index;
			index++;
			p = q;
		}
	}

	map<int, int> m;
	for(int i = 0; i < v.size(); i++)
	{
		assert(v[i] != -1);
		int x = im[i];
		m.insert(PI(x, v[i]));
	}

	return m;
}

int shseg::draw(const string & file, int xsize)
{
	ofstream fout(file.c_str());

	map<int, int> m1 = sort_vertices(xsize, true);
	map<int, int> m2 = sort_vertices(xsize, false);

	double px = 0;
	double py = 0;

	if(m1.size() < m2.size()) px = (m2.size() - m1.size()) / 2;
	else py = (m1.size() - m2.size()) / 2;
	
	draw_header(fout);

	fout<<"\\def\\glen{1.5cm}\n";
	fout<<"\\def\\ylen{6cm}\n";

	map<int, int>::iterator it;
	// vertices
	char sx[1024];
	char sy[1024];
	for(it = m1.begin(); it != m1.end(); it++)
	{
		int x = it->first;
		int a = ie[x];
		int p = it->second;
		sprintf(sx, "s%d", x);
		string cx = "colx";
		if(sa.find(a) != sa.end()) cx = "colc";
		if(sb.find(a) != sb.end()) cx = "colb";

		fout<<"\\node[mycircle, \\"<<cx.c_str()<<", draw, label = above:"<<a<<"] ("<<sx<<") at ("<<p + px<<" * \\glen, 0 * \\ylen) {"<<x<<"};\n";
	}

	for(it = m2.begin(); it != m2.end(); it++)
	{
		int x = it->first;
		int a = ie[x];
		int p = it->second;
		sprintf(sx, "s%d", x);
		string cx = "colx";
		if(sa.find(a) != sa.end()) cx = "colc";
		if(sb.find(a) != sb.end()) cx = "colb";

		fout<<"\\node[mycircle, \\"<<cx.c_str()<<", draw, label = below:"<<a<<"] ("<<sx<<") at ("<<p + py<<" * \\glen, -1 * \\ylen) {"<<x<<"};\n";
	}

	// edges in eb
	for(set<PI>::iterator it = eb.begin(); it != eb.end(); it++)
	{
		int a = it->first < it->second ? it->first : it->second;
		int b = it->first > it->second ? it->first : it->second;
		int x = ei[a];
		int y = ei[b];
		sprintf(sx, "s%d", x);
		sprintf(sy, "s%d", y);
		if(a < xsize && b >= xsize)
		{
			fout<<"\\draw[\\colb] ("<<sx<<") -- ("<<sy<<");\n";
		}
		else if(a < xsize && b < xsize)
		{
			fout<<"\\draw[\\colb, bend left = 30] ("<<sx<<") to ("<<sy<<");\n";
		}
		else if(a >= xsize && b >= xsize)
		{
			fout<<"\\draw[\\colb, bend right = 30] ("<<sx<<") to ("<<sy<<");\n";
		}
		else assert(false);
	}

	// edges in ea
	for(set<PI>::iterator it = ea.begin(); it != ea.end(); it++)
	{
		int a = it->first < it->second ? it->first : it->second;
		int b = it->first > it->second ? it->first : it->second;
		int x = ei[a];
		int y = ei[b];
		sprintf(sx, "s%d", x);
		sprintf(sy, "s%d", y);
		string c = "colc";
		string s = "";
		if(contain(a / 2, b / 2)) s = "thick";
		fout<<"\\draw["<<s.c_str()<<",\\"<<c.c_str()<<"] ("<<sx<<") -- ("<<sy<<");\n";
	}
		
	set<int> sdup;
	for(int i = 0; i < vdup.size(); i++)
	{
		for(int j = 0; j < vdup[i].size(); j++)
		{
			int x = vdup[i][j];
			if(sdup.find(x) != sdup.end()) continue;
			sdup.insert(x);
		}
	}

	// edges for genes
	for(int i = 0; i < ge.size(); i++)
	{
		int x = ge[i].first;
		int y = ge[i].second;
		sprintf(sx, "s%d", x);
		sprintf(sy, "s%d", y);
		string c = "cold";
		if(sdup.find(i) != sdup.end()) c = "colf";
		string p = "above";
		if(ie[x] < xsize) fout<<"\\draw[line width = 0.07cm, \\"<<c.c_str()<<", ->] ("<<sx<<") to node[label=above:"<<gf[i]<<"]{} ("<<sy<<");\n";
		else fout<<"\\draw[line width = 0.07cm, \\"<<c.c_str()<<", ->] ("<<sx<<") to node[label=below:"<<gf[i]<<"]{} ("<<sy<<");\n";
	}

	// edges for adjacencies
	for(int i = 0; i < adj.size(); i++)
	{
		int j = adj[i];
		if(j == -1) continue;
		if(i > j) continue;
		sprintf(sx, "s%d", i);
		sprintf(sy, "s%d", j);
		int a = ie[i];
		int b = ie[j];
		if(a < xsize && b >= xsize)
		{
			fout<<"\\draw[very thick, \\cola] ("<<sx<<") -- ("<<sy<<");\n";
		}
		else if(a < xsize && b < xsize)
		{
			if((int)fabs(m1[i] - m1[j]) != 1) fout<<"\\draw[very thick, \\cola, bend right = 30] ("<<sx<<") to ("<<sy<<");\n";
			else fout<<"\\draw[very thick, \\cola] ("<<sx<<") -- ("<<sy<<");\n";
		}
		else if(a >= xsize && b >= xsize)
		{
			if((int)fabs(m2[i] - m2[j]) != 1) fout<<"\\draw[very thick, \\cola, bend left = 30] ("<<sx<<") to ("<<sy<<");\n";
			else fout<<"\\draw[very thick, \\cola] ("<<sx<<") -- ("<<sy<<");\n";
		}
		else assert(false);

		//fout<<"\\draw[very thick, \\cola, bend left = 30] ("<<sx<<") to ("<<sy<<");\n";
	}

	draw_footer(fout);
	fout.close();
	return 0;
}
