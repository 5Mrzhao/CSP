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
void feGRASS(int *ai_in, int *aj_in, double *av_in,
			int N_in, int m1, int m2, int m3, int m4, int* insparsifier_in, double alpha_in);
typedef struct matGraph {
    int n;
    double time_feGRASS;
    std::vector<edge> edges;
    matGraph(int n, int m) : n(n) {edges.reserve(m);}
    void sparsify1(double p, int m1, int m2, int m3, int m4,
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
        feGRASS(u2, v2, w2, n, m1, m2, m3, m4, in2, p);
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
        for(int tar = m1 + m2; i < tar; i++)
        {
            if(in2[i]!=0)
            {
                u3.push_back(edges[i].u);
                v3.push_back(edges[i].v);
                w3.push_back(edges[i].w);
            }
        }
        for(int i1, j1, tar = m1 + m2 + m3; i < m1 + m2 + m3; i++, tar++){
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
}matGraph;