#include <iostream>
#include <vector>
#include <ctime>
#include <math.h>
#include <algorithm>
#include <string.h>
#include <stack>
#include <unordered_map>
#include <map>
#include <sys/time.h>
#include "./structs.h"

#define VALUE_TYPE double
using namespace std;
#define VALUE_TYPE double
using namespace std;
// void feGRASS2(int *ai_in, int *aj_in, double *av_in,
// 			int N_in, int m1, int m2, int m3, int* insparsifier_in, double alpha_in);
typedef struct matGraph2{
    int n;
    double time_feGRASS;
    std::vector<edge> edges;
    matGraph2(int n, int m) : n(n) {edges.reserve(m);}
    void sparsify1(double p, int m1, int m2, int m3,
                std::vector<int>& Row_Num_end_fegrass, std::vector<int>& Col_Num_end_fegrass, std::vector<VALUE_TYPE>& Val_Num_end_fegrass)
    {   
        int m = edges.size();
        int* u2 = new int[m];
        int* v2 = new int[m];
        double* w2 = new double[m];
        int* in2 = new int[m]; memset(in2, 0, sizeof(int) * m);
        for (int i = 0; i < m; i++)
        {
            u2[i] = edges[i].u;
            v2[i] = edges[i].v;
            w2[i] = fabs(edges[i].w);
        }
        feGRASS2(u2, v2, w2, n, m1, m2, m3, in2, p);
        std::vector<int> u3, v3;
        std::vector<VALUE_TYPE> w3;
        int i, all;
        m = m1 + m2 + m3;
        u3.reserve(m), v3.reserve(m), w3.reserve(m);
        for(i = 0; i < m1; i++)
        {
            if(in2[i]!=0)
            {
                u3.push_back(edges[i].u);
                v3.push_back(edges[i].v);
                w3.push_back(edges[i].w);
                u3.push_back(edges[i].v);
                v3.push_back(edges[i].u);
                w3.push_back(edges[i].w);
            }
        }
        for(int i1, j1, tar = m1 + m2; i < m1 + m2; i++, tar++)
        {
            if(in2[i]!=0)
            {
                i1 = edges[i].u; j1 = edges[i].v;
                u3.push_back(edges[i].u);
                v3.push_back(edges[i].v);
                w3.push_back(edges[i].w);
                u3.push_back(edges[tar].u);
                v3.push_back(edges[tar].v);
                w3.push_back(edges[tar].w);
            }
        }
        all = u3.size();
        Row_Num_end_fegrass.reserve(all + n);
        Col_Num_end_fegrass.reserve(all + n);
        Val_Num_end_fegrass.reserve(all + n);
        for(int i = 0; i < all; i++)
        {
            Row_Num_end_fegrass.push_back(u3[i]);
            Col_Num_end_fegrass.push_back(v3[i]);
            Val_Num_end_fegrass.push_back(w3[i]);
        }
        delete[] u2;
        delete[] v2;
        delete[] w2;
        delete[] in2;
        u3.clear();
        v3.clear();
        w3.clear();
    }
    void addEdge(int u, int v, VALUE_TYPE w)
    {
        edges.push_back(edge(u, v, w));
    }
    Graph gra;
    bfsnode* bfsque;

    int* ai, * aj, * edgeperm, * deg, * layer, * LCA, * insparsifier, * fa, * vis, * ans, * parent;
    int SIZE, M, M1, M2, M3, maxnode = 1, totedges = 0, treeedges;
    double* av, * vol, * dist, * score, * effwt;
    double totstretch, avstretch, scale = -1.0;

    inline int getfa1(int x) { return fa[x] == x ? x : fa[x] = getfa1(fa[x]); }
    void BFS1()
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
    void Kruscal1()
    {

        double maxdeg, sumlayer = 1, scale = 1; int i, j;
        BFS1();
        for (int l = 0; l < M; l++)
        {			
            i = ai[l]; j = aj[l];
            maxdeg = (double)(deg[i] > deg[j] ? deg[i] : deg[j]);
            sumlayer = (double)(layer[i] + layer[j]);
            scale = log (maxdeg) /  (sumlayer);
            effwt[l] = av[l] * scale;
        }
        sort(edgeperm, edgeperm + M, [this](int i,int j){return this->effwt[i] > this->effwt[j];});
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
        treeedges = count; totedges = treeedges;
        dis.destroy();
    }
    void BFSmark1(int* mark, int i, int layer, int* visited, int rank, int nodenum)
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

    void addedges1(double fac, int sigma)
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
        sort(edgeperm, edgeperm + (ptr1), [this](int i,int j){return this->effwt[i] > this->effwt[j];});
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
                BFSmark1(deg, i, sigma, layer, count, nodenum);
                BFSmark1(deg, j, sigma, layer, count, nodenum);
                insparsifier[ind] = true;
            }
            ptr++;
        }
        totedges += count;
        //cout << totedges << " edges recovered" << endl;
        return;
    }
    void dfs1(int start, Graph& Tree)
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
                hasson = true;
                len = 1.0 / (av[ind]);
                dist[child] = dist[node] + len;
                deg[node]++;
                dfsstack.push(dfsnode{ child, node });
                break;
            }
            if (!hasson)
            {
                dfsstack.pop();
                for (int k = gra.ptr[node]; k < gra.ptr[node + 1]; k++)
                {
                    if (layer[gra.edgelists[k].node])
                    {
                        int ind = gra.edgelists[k].ind;
                        LCA[ind] = getfa1(gra.edgelists[k].node);
                    }
                }
                fa[node] = from;
                layer[node] = 1;
            }
        }
    }
    void tarjan1(int root, Graph& Tree)
    {
        totstretch = 0;
        for (int i = 0; i < SIZE; i++)  fa[i] = i;
        memset(layer, 0, sizeof(int) * (SIZE));
        for (int i = 0; i < SIZE; i++)
        {
            if (layer[i] == 0)
            {
                dist[i] = 0;
                dfs1(i, Tree);
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
    void calculateresistance1()
    {
        tarjan1(0, gra);
    }
    void initialize1() {
        deg = (int*)malloc(sizeof(int) * (SIZE)); memset(deg, 0, sizeof(int) * (SIZE));
        layer = (int*)malloc(sizeof(int) * (SIZE)); memset(layer, -1, sizeof(int) * (SIZE));
        dist = (double*)malloc(sizeof(double) * (SIZE)); memset(dist, 0, sizeof(double) * (SIZE));
        bfsque = (bfsnode*)malloc(sizeof(bfsnode) * (SIZE));
        fa = (int*)malloc(sizeof(int) * (SIZE)); memset(fa, 0, sizeof(int) * (SIZE));

        edgeperm = (int*)malloc(sizeof(int) * (M)); LCA = (int*)malloc(sizeof(int) * (M));
        effwt = (double*)malloc(sizeof(double) * (M));
        memset(insparsifier, 0, sizeof(int) * (M));

        for (int l = 0; l < M; l++) { edgeperm[l] = l; }
    }
    void constructgraph1()
    {
        gra = Graph(SIZE, M);
        int i, j;
        for (int l = 0; l < M; l++)
        {
            i = ai[l]; j = aj[l];
            deg[i]++; deg[j]++;
        }
        gra.edgelists = (Edge*)malloc(sizeof(Edge) * (2 * M));
        gra.ptr = (int*)malloc(sizeof(int) * (SIZE + 1));
        int* ptrs = (int*)malloc(sizeof(int) * (SIZE));
        maxnode = 0; gra.ptr[0] = 0; ptrs[0] = 0;
        int sum = 0;
        for (int l = 1; l < SIZE; l++)
        {
            if (deg[l] > deg[maxnode]) maxnode = l;
            sum += deg[l - 1];
            gra.ptr[l] = sum;
            ptrs[l] = sum;
        }
        sum += deg[SIZE - 1];
        gra.ptr[SIZE] = sum;
        for (int l = 0; l < M; l++)
        {
            i = ai[l]; j = aj[l];
            gra.edgelists[ptrs[i]] = Edge{ j,l };
            gra.edgelists[ptrs[j]] = Edge{ i,l };
            ptrs[i]++; ptrs[j]++;
        }
        free(ptrs);
    }

    void freeMemory1()
    {
        gra.destroy();
        free(edgeperm); free(LCA);
        free(deg); free(layer); free(bfsque); free(dist);
    }
    void Change1()
    {
        for(int a = M1, p = M1 + M2; a < M; ++a, ++p)
        {
            if(av[p] < av[a]) av[a] = av[p];
        }
    }
    void feGRASS2(int *ai_in, int *aj_in, double *av_in, int N_in, int m1, int m2, int m3, int* insparsifier_in, double alpha_in)
    {
        ai = ai_in; aj = aj_in; av = av_in; M = m1 + m2; SIZE = N_in; insparsifier = insparsifier_in; double fac = alpha_in;
        M1 = m1, M2 = m2, M3 = m3;
        initialize1();
        constructgraph1();
        Change1();
        Kruscal1();
        calculateresistance1();
        addedges1(fac, 1);
        freeMemory1();
        return;
    }
}matGraph2;