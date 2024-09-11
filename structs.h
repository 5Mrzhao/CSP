#pragma once
#include <cstdlib>
struct Edge
{
	int node;
	int ind;
};
class Graph
{
public:
	int n; int m;
	Edge* edgelists;
	int* ptr;
	Graph(int n, int m) :n(n), m(m) { }
	Graph() { m = 0; }
	void destroy() { free(edgelists); free(ptr); }
};
class DisjSet
{
public:
	int* parent;
	int* rank;
	DisjSet(int max_size)
	{
		parent = (int*)malloc(sizeof(int) * max_size);
		rank = (int*)malloc(sizeof(int) * max_size);
		for (int i = 0; i < max_size; i++) parent[i] = i;
	}
	int find(int x)
	{
		int p = parent[x];
		if (parent[p] != p)
		{
			p = find(p);
			parent[x] = p;
		}
		return p;
	}
	void to_union(int x1, int x2)
	{
		int f1 = find(x1);
		int f2 = find(x2);
		if (f1 == f2) return;
		if (rank[f1] < rank[f2]) parent[f1] = f2;
		else
		{
			parent[f2] = f1;
			if (rank[f1] == rank[f2]) ++rank[f1];
		}
	}
	bool is_same(int e1, int e2) { return find(e1) == find(e2); }
	void destroy() { free(parent); free(rank); }
};
struct dfsnode { int node, from; };
struct bfsnode { int node, layer; };
typedef struct edge
{
    int u, v;
    double w;
    edge(int u, int v, double w) : u(u), v(v), w(w){}
}edge;