#include "mygraph.h"

int calc_shortest_cycle(int start, vector<int> & list, const dgraph & gr, const vector<bool> & del)
{
	assert(start >= 0 && start < num_vertices(gr));

	assert(del[start] == false);

	// queque, stored those that are not expanded
	vector<int> q;

	// store whether this node is touched
	vector<bool> v;
	v.assign(num_vertices(gr), false);

	// store whether this node is searched, its level and parent if searched
	vector<int> level;
	vector<int> parent;

	level.assign(num_vertices(gr), -1);
	parent.assign(num_vertices(gr), -1);

	q.push_back(start);
	v.at(start) = true;
	level.at(start) = 0;
	parent.at(start) = -1;

	int qi = 0;

	bool found = false;

	while(qi < q.size())
	{
		int index = q.at(qi);

		assert(v.at(index)==true);
		assert(del[index] == false);

		//MIV::iterator it;
		//for(it = nodes.at(index)->out.begin(); it != nodes.at(index)->out.end(); it++)
		out_edge_iterator_d ei1, ei2;
		for(tie(ei1, ei2) = out_edges(index, gr); ei1 != ei2; ei1++)
		{
			//assert(del[it->first] == false);
			int s = source(*ei1, gr);
			int t = target(*ei1, gr);
			assert(s == index);
			if(del[t] == true) continue;

			if(t == start)
			{
				list.clear();
				int p = index;
				list.push_back(p);

				while(true)
				{
					p = q.at( parent.at(p) );
					list.push_back(p);
					if(p == start) break;
				}

				found = true;
				break;
			}

			if(v.at(t) == true) continue;

			v.at(t) = true;
			level.at(t) = level.at(index) + 1;
			parent.at(t) = qi;
			q.push_back(t);
		}

		if(found == true) break;

		qi++;
	}
	return 0;
}

int calc_cycles(vector< vector<int> > & cycles, const dgraph & gr)
{
	vector<bool> del;
	del.assign(num_vertices(gr), false);

	cycles.clear();
	for(int i = 0; i < num_vertices(gr); i++)
	{
		if(del[i] == true) continue;
		vector<int> v;
		calc_shortest_cycle(i, v, gr, del);
		if(v.size() > 0)
		{
			cycles.push_back(v);
			del[i] = true;
		}
	}

	//printf("found %5d cycles\n", (int)(cycles.size()));
	return 0;
}
