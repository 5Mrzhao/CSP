#include "./Sfegrass.hpp"
using namespace std;

Graph gra;
bfsnode* bfsque;

int *ai, *aj, *edgeperm, *deg, *layer, *LCA, *insparsifier, *fa, *vis, *ans, *parent;
int SIZE, M, M1, M2, M3, M4, maxnode = 1, totedges = 0, treeedges;
double* av, * dist, * score, * effwt, *vol;
double totstretch, avstretch, scale = -1.0;

inline int getfa(int x) { return fa[x] == x ? x : fa[x] = getfa(fa[x]); }
bool compareByEffWt(int i, int j){return effwt[i] > effwt[j];}
bool compareByWt(Edge i, Edge j) { return av[i.ind] > av[j.ind]; }
bool compareByDeg(int i, int j) { return deg[i] > deg[j]; }

void initialize() {
	deg = (int*)malloc(sizeof(int) * (SIZE)); memset(deg, 0, sizeof(int) * (SIZE));
	vol = (double *)malloc(sizeof(double) * (SIZE)); memset(vol, 0, sizeof(double) * SIZE);

	layer = (int*)malloc(sizeof(int) * (SIZE)); memset(layer, -1, sizeof(int) * (SIZE));
	dist = (double*)malloc(sizeof(double) * (SIZE)); memset(dist, 0, sizeof(double) * (SIZE));
	bfsque = (bfsnode*)malloc(sizeof(bfsnode) * (SIZE));
	fa = (int*)malloc(sizeof(int) * (SIZE)); memset(fa, 0, sizeof(int) * (SIZE));

	edgeperm = (int*)malloc(sizeof(int) * (M)); LCA = (int*)malloc(sizeof(int) * (M));
	effwt = (double*)malloc(sizeof(double) * (M));
	memset(insparsifier, 0, sizeof(int) * (M));

	for(int l = 0; l < M; l++) { edgeperm[l] = l; }
}
void constructgraph()
{
	gra = Graph(SIZE, M);
	int i, j, l;
	double val;
	for (l = 0; l < M1; l++)
	{
		i = ai[l], j = aj[l];
		++deg[i], ++deg[j];
		vol[i] += av[l];
		vol[j] += av[l];
	}
	for (; l < M1 + M2; l++)
	{
		i = ai[l], j = aj[l];
		++deg[i], ++deg[j];
		vol[i] += av[l] / 2.0;
	}
	for (int p = M; l < M; l++, p++)
	{
		i = ai[l], j = aj[l];
		++deg[i], ++deg[j];
		val = (av[l] + av[p]) / 2.0;
		vol[i] += val;
		vol[j] += val;
	}
	gra.edgelists = (Edge*)malloc(sizeof(Edge) * (2 * M));
	gra.ptr = (int*)malloc(sizeof(int) * (SIZE + 1));
	int* ptrs = (int*)malloc(sizeof(int) * (SIZE));
	maxnode = 0; gra.ptr[0] = 0; ptrs[0] = 0;
	int sum = 0;
	for (l = 1; l < SIZE; l++)
	{
		if (deg[l] > deg[maxnode]) maxnode = l;
		sum += deg[l - 1];
		gra.ptr[l] = sum;
		ptrs[l] = sum;
	}
	sum += deg[SIZE - 1];
	gra.ptr[SIZE] = sum;
	for (l = 0; l < M; l++)
	{
		i = ai[l]; j = aj[l];
		gra.edgelists[ptrs[i]] = Edge{ j,l };
		gra.edgelists[ptrs[j]] = Edge{ i,l };
		ptrs[i]++;
		ptrs[j]++;
	}
	free(ptrs);
}
void freeMemory()
{
	gra.destroy();
	free(edgeperm); free(LCA);
	free(deg); free(layer); free(bfsque); free(dist);
}

void BFS()
{
	int head = 0, tail = 0, node, curlayer, i = maxnode, comp = 0;
	bfsque[tail] = bfsnode{ i,0 }; tail++;
	layer[i] = 0;
	while (head < tail)
	{
		node = bfsque[head].node; curlayer = bfsque[head].layer; head++;
		for (int k = gra.ptr[node]; k < gra.ptr[node + 1]; k++)
		{
			if (layer[gra.edgelists[k].node] < 0)
			{
				bfsque[tail++] = bfsnode{ gra.edgelists[k].node,curlayer + 1 };
				layer[gra.edgelists[k].node] = curlayer + 1;
			}
		}
	}
	for (int i = 0; i < SIZE; i++)
	{
		if (layer[i] >= 0) continue;
		head = 0; tail = 0;
		bfsque[tail] = bfsnode{ i,0 }; tail++;
		layer[i] = 0;
		while (head < tail)
		{
			node = bfsque[head].node; curlayer = bfsque[head].layer; head++;
			for (int k = gra.ptr[node]; k < gra.ptr[node + 1]; k++)
			{
				if (layer[gra.edgelists[k].node] < 0)
				{
					bfsque[tail++] = bfsnode{ gra.edgelists[k].node,curlayer + 1 };
					layer[gra.edgelists[k].node] = curlayer + 1;
				}
			}
		}
	}
}
void Kruscal()
{
	double sumw, sumv, sumlayer = 1, scale = 1;
	int i, j, s = 0, e1 = M1 + M2, e2 = M;
	BFS();
	for (; s < e1; s++)
	{
		i = ai[s], j = aj[s];
		sumw = av[s] * 2;
		sumv = vol[i] + vol[j];
		sumlayer = (double)(layer[i] + layer[j]) * 2;
		scale = log(sumv) / (sumlayer);
		effwt[s] = sumw * scale;
	}
	for (int p = M; s < M; s++, p++)
	{
		i = ai[s], j = aj[s];
		sumw = av[s] + av[p];
		sumv = vol[i] + vol[j];
		sumlayer = (double)(layer[i] + layer[j]) * 2;
		scale = log(sumv) / (sumlayer);
		effwt[s] = sumw * scale;
	}
	sort(edgeperm, edgeperm + M, compareByEffWt);
	DisjSet dis(SIZE);
	int ind;
	int count = 0;
	for (int l = 0; l < M; l++)
	{
		if (count >= SIZE - 1) { break; }
		ind = edgeperm[l];
		i = ai[ind]; j = aj[ind];
		if (!dis.is_same(i, j))
		{
			insparsifier[ind] = true;
			count++;
			dis.to_union(i, j);
		}
	}
	//printf("%d,\n",count);
	treeedges = count; totedges = treeedges;
	dis.destroy();
	//cout << "Kruscal end, " << count << " edges in spanning tree, " << SIZE - count << " connected components." << endl;
}
void Change()
{
	int s = M1 + M2, e = M, i, j;
	for(int p = M; s < e; ++s, ++p)
	{
		if(insparsifier[s]) continue;
		if(av[s] < av[p]) av[s] = av[p];
	}
}

void dfs(int start, Graph& Tree)
{
	stack<dfsnode> dfsstack;
	dfsstack.push(dfsnode{ start, -1 });
	int ptr, count = 0;
	memset(deg, 0, sizeof(int) * (SIZE));
	int node, child, from; int ind; double len;
	bool hasson;
	while (!dfsstack.empty())
	{
		node = dfsstack.top().node;
		from = dfsstack.top().from;
		hasson = false;
		int length = Tree.ptr[node + 1] - Tree.ptr[node];
		int start = Tree.ptr[node];
		for (; deg[node] < length;)
		{
			ptr = deg[node];
			child = Tree.edgelists[start + ptr].node;
			ind = Tree.edgelists[start + ptr].ind;
			if (deg[child] != 0 || !insparsifier[ind])
			{
				deg[node]++;
				continue;
			}
			//cout << "OK" << endl;
			hasson = true;
			len = 1.0 / (av[ind]);
			dist[child] = dist[node] + len;
			deg[node]++;
			dfsstack.push(dfsnode{ child, node });
			break;
		}
		//cout << "OK1" << endl;
		if (!hasson)
		{
			dfsstack.pop();
			for (int k = gra.ptr[node]; k < gra.ptr[node + 1]; k++)
			{
				if (layer[gra.edgelists[k].node])
				{
					int ind = gra.edgelists[k].ind;
					LCA[ind] = getfa(gra.edgelists[k].node);
				}
			}
			fa[node] = from;
			layer[node] = 1;
		}
		//cout << "OK2" << endl;
	}
}
void tarjan(int root, Graph& Tree)
{
	totstretch = 0;
	for (int i = 0; i < SIZE; i++)  fa[i] = i;
	memset(layer, 0, sizeof(int) * (SIZE));
	for (int i = 0; i < SIZE; i++)
	{
		if (layer[i] == 0)
		{
			dist[i] = 0;
			dfs(i, Tree);
		}
	}
	int i, j;
	double w, R, longest = 0;
	totstretch = 0;
	for (int l = 0; l < M; l++)
	{
		i = ai[l];
		j = aj[l];
		w = av[l];
		R = (dist[i] + dist[j] - 2 * dist[LCA[l]]);
		effwt[l] = w * R;
		totstretch += effwt[l];
	}
	avstretch = totstretch / M;
	free(fa);
}
void calculateresistance()
{
	tarjan(0, gra);
}

void BFSmark(int* mark, int i, int layer, int* visited, int rank, int nodenum)
{
	int head = 0, tail = 0, node, flag = 0, curlayer;
	bfsque[tail] = bfsnode{ i,0 }; tail++;
	while (head < tail && tail < nodenum)
	{
		node = bfsque[head].node; curlayer = bfsque[head].layer; head++;
		if (curlayer >= layer) { break; }
		for (int k = gra.ptr[node]; k < gra.ptr[node + 1]; k++)
		{
			if (!insparsifier[gra.edgelists[k].ind] || visited[gra.edgelists[k].node] >= rank) continue;
			visited[gra.edgelists[k].node] = rank;
			bfsque[tail] = bfsnode{ gra.edgelists[k].node,curlayer + 1 };
			tail++;
		}
	}
	for (int j = 0; j < tail; j++) { mark[bfsque[j].node]++; }
}
void addedges(double fac, int sigma)
{
	int ptr1 = 0, ptr2 = M - 1;
	for (int l = 0; l < M; l++)
	{
		if (insparsifier[l]) {
			edgeperm[ptr2] = l; ptr2--;
		}
		else {
			edgeperm[ptr1] = l; ptr1++;
		}
	}
	sort(edgeperm, edgeperm + (ptr1), compareByEffWt);
	free(effwt);
	int uplink = 1;
	int num = int(fac * double(SIZE));
	int nodenum = (int)(0.6 / fac);
	int count = 0;
	int ptr = 0, ind, bdry = M / 2;
	// printf("%d\n",M);
	memset(deg, 0, sizeof(int) * (SIZE));
	memset(layer, 0, sizeof(int) * (SIZE));
	int i, j;

	while (1) {
		// printf("%d ?? %d\n",count, num);
		if (count >= num) break;
		// printf("%d ? %d\n",ptr, bdry);
		if (ptr >= bdry) {
			ptr = 0; uplink++;
			if (uplink > SIZE) { break; }
		}
		ind = edgeperm[ptr];
		i = ai[ind]; j = aj[ind];
		if (!insparsifier[ind] && (deg[i] < uplink && deg[j] < uplink))
		{
			count++;
			BFSmark(deg, i, sigma, layer, count, nodenum);
			BFSmark(deg, j, sigma, layer, count, nodenum);
			insparsifier[ind] = true;
		}
		ptr++;
	}
	totedges += count;
	//cout << totedges << " edges recovered" << endl;
	return;
}

void feGRASS(int *ai_in, int *aj_in, double *av_in, 
			int N_in, int M1_in, int M2_in, int M3_in, int M4_in, int* insparsifier_in, double alpha_in)
{
	ai = ai_in; aj = aj_in; av = av_in;
	M1 = M1_in, M2 = M2_in, M3 = M3_in, M4 = M4_in;
	M = M1 + M2 + M3; SIZE = N_in;
	insparsifier = insparsifier_in;
	initialize();
	constructgraph();
	Kruscal();
	Change();
	calculateresistance();
	addedges(alpha_in, 1);
	freeMemory();
	return;
}