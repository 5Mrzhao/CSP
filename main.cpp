#include <bits/stdc++.h>
#include "./Sfegrass.hpp"
#include "./fegrass_1030.hpp"
#include "./HB_ADD.hpp"
using namespace std;

// string saveMtxPath1 = "/home/zyx/getData/saveMtx/";
// string saveMtxPath2 = "/home/zyx/getData/saveMtx/";
// string saveMtxPath3 = "/home/zyx/getData/saveMtx/";
// string saveMtxPath4 = "/home/zyx/getData/saveMtx/name.csv";
string saveMtxPath1 = "/home/zyx/ZYX_iccad24/getData/saveMtx/";
string saveMtxPath2 = "/home/zyx/ZYX_iccad24/getData/saveMtx/";
string saveMtxPath3 = "/home/zyx/ZYX_iccad24/getData/saveMtx/";
string saveMtxPath4 = "/home/zyx/ZYX_iccad24/getData/saveMtx/name.csv";
void getMtxName(char str[])
{
    int tar = 0, p = 0;
    while(str[tar] != '.') ++tar;
    while(str[tar] != '/') --tar;
    ++tar;
    while(1)
    {
        str[p] = str[tar];
        if(str[tar] == '.') break;
        ++p, ++tar;
    }
    str[p] = '\0'; 
    saveMtxPath1 += str; saveMtxPath1 += "";
    saveMtxPath2 += str; saveMtxPath2 += "_cspa";
    saveMtxPath3 += str; saveMtxPath3 += "_gpscpa";
    return;
}
void getMtxName1(char str[], char name[])
{
    int tar = 0, p = 0;
    while(str[tar] != '.') ++tar;
    while(str[tar] != '/') --tar;
    ++tar;
    while(1)
    {
        str[p] = str[tar];
        name[p] = str[tar];
        if(str[tar] == '.') break;
        ++p, ++tar;
    }
    str[p] = '\0'; 
    name[p] = '\0';
    name = str;
    saveMtxPath1 += str; saveMtxPath1 += "";
    saveMtxPath2 += str; saveMtxPath2 += "_csp";
    saveMtxPath3 += str; saveMtxPath3 += "_gpscp";
    return;
}
void saveMtx(string filename, int &n, std::vector<int>& Row_Num, std::vector<int>& Col_Num, std::vector<VALUE_TYPE>& Val_Num)
{
    ofstream fout(filename, std::ios::out);
    // fout<<"%%MatrixMarket matrix coordinate real general"<<endl;
    fout << n << ' ' << n << ' ' << Val_Num.size() << endl;
    for(int a = 0; a < Val_Num.size(); ++a) fout << Row_Num[a]+1 << ' ' << Col_Num[a]+1 << ' ' << setprecision(15) << Val_Num[a] << endl;
    fout.close();
}
void saveMtx1(string filepath, char name[])
{
    // ofstream fout(filepath, std::ios::out);
    // fout << name << endl;
    // fout.close();
    // string filepath = "/home/zyx/getData/saveMtx/name.csv";
    freopen(filepath.c_str(), "a", stdout);
    printf("%s\n",name);
    fclose(stdout);
}
void readMtx(char *filename, int &n, int &m, std::vector<int>& Row_Num, std::vector<int>& Col_Num, std::vector<VALUE_TYPE>& Val_Num)
{
    ifstream fin(filename, std::ios::in);
    while(fin.peek() == '%') fin.ignore(2048, '\n');
    if(!fin.is_open())
    {
        cout << "ERROR filename" << endl;
        return;
    }
    int i, j, m1=0;
    double val;
    fin >> n >> i >> m;
    Row_Num.resize(m), Col_Num.resize(m), Val_Num.resize(m);
    cout << filename << endl;
    cout << n << ' ' << n << ' ' << m << endl;
    if(n != i)
    {
        cout << "ERROR n1!=n2" << endl;
        fin.close();
        return;
    }
    for(int a = 0; a < m; ++a)
    {
        fin >> i >> j >> val;
        --i;
        --j;
        Row_Num[a] = i, Col_Num[a] = j, Val_Num[a] = val;
    }
    fin.close();
    return;
}

/*
0: 不对称边
1: 完全对称边
>=2: 值不对称边
*/
void getIssystem(int n, int m, std::vector<int>& Row_Num, std::vector<int>& Col_Num, std::vector<VALUE_TYPE>& Val_Num,
                    std::vector<int> &Issystem)
{
    int i, j, tar;
    VALUE_TYPE w;
    vector<unordered_map<int, int>> vecm; vecm.resize(n);
    for(int l = 0; l < m; ++l)
    {
        i = Row_Num[l], j = Col_Num[l];
        (vecm[i])[j] = l + 1;
    }
    for(int l = 0; l < m; ++l)
    {
        if(Issystem[l] != 0) continue;
        i = Row_Num[l], j = Col_Num[l];
        if((vecm[j])[i] != 0)
        {
            tar = (vecm[j])[i] - 1;
            if(Val_Num[l] != Val_Num[tar])
            {
                if(abs(Val_Num[l]) > abs(Val_Num[tar]))
                {
                    Issystem[l] = 2+l;
                    Issystem[tar] = 2+l;
                }
                else
                {
                    Issystem[l] = 2+tar;
                    Issystem[tar] = 2+l;
                }
            }
            else
            {
                Issystem[l] = 1;
                Issystem[tar] = 1;
            }
        }
    }
}

/*
m1: 完全对称边下三角
m2: 非对称边
m3: 值不对称边下三角
m4: 值不对称边上三角
*/
void My_Matrix_Sys1(int n, int nnzR, int &m1, int &m2, int &m3, int &m4,
                    std::vector<int>& Row_Number_coo, std::vector<int>& Col_Number_coo, std::vector<VALUE_TYPE>& Val_Number_coo,              
                    std::vector<int>& Row_Number_my_sys_for_fegrass_angle, std::vector<int>& Col_Number_my_sys_for_fegrass_angle, std::vector<VALUE_TYPE>& Val_Number_my_sys_for_fegrass_angle,
                    std::vector<int>& Issystem)
{
    struct timeval t1, t2;
    int i,j,k, p1, p2, p3, p4, tar, tar1, temp; m1 = m2 = m3 = m4 = 0;
    VALUE_TYPE w;
    double time_MY_SYS = 0;
    gettimeofday(&t1, NULL);
    for(i = 0; i < nnzR; ++i)
    {
        if(Issystem[i] == 0 && Row_Number_coo[i] != Col_Number_coo[i]) ++m2;
        else if(Issystem[i] == 1 && Row_Number_coo[i] > Col_Number_coo[i]) ++m1;
        else if(Issystem[i] >= 2 && Row_Number_coo[i] > Col_Number_coo[i]) ++m3, ++m4;
    }
    Row_Number_my_sys_for_fegrass_angle.resize(m1+m2+m3+m4);
    Col_Number_my_sys_for_fegrass_angle.resize(m1+m2+m3+m4);
    Val_Number_my_sys_for_fegrass_angle.resize(m1+m2+m3+m4);
    p1 = 0, p2 = m1, p3 = m1 + m2, p4 = m1 + m2 + m3;
    for(int a = 0; a < nnzR; ++a)
    {
        tar = -1;
        if(Row_Number_coo[a] == Col_Number_coo[a]) continue;
        i = Row_Number_coo[a], j = Col_Number_coo[a], w = Val_Number_coo[a];
        if(Issystem[a] == 0 && Row_Number_coo[a] != Col_Number_coo[a]) tar = p2++;
        else if(Issystem[a] == 1 && Row_Number_coo[a] > Col_Number_coo[a]) tar = p1++;
        else if(Issystem[a] >= 2 && Row_Number_coo[a] > Col_Number_coo[a])
        {
            tar = p3++;
            tar1 = p4++;
            Row_Number_my_sys_for_fegrass_angle[tar1] = j;
            Col_Number_my_sys_for_fegrass_angle[tar1] = i;
            Val_Number_my_sys_for_fegrass_angle[tar1] = Val_Number_coo[Issystem[a]-2];
        }
        else continue;
        Row_Number_my_sys_for_fegrass_angle[tar] = i;
        Col_Number_my_sys_for_fegrass_angle[tar] = j;
        Val_Number_my_sys_for_fegrass_angle[tar] = w;
    }
    gettimeofday(&t2, NULL);
    time_MY_SYS = (t2.tv_sec - t1.tv_sec) * 1000.0 + (t2.tv_usec - t1.tv_usec) / 1000.0;
    // printf("time_MY_SYS_for_fegrass=%lf ms\n", time_MY_SYS);
}

/*
m1: 完全对称边下三角
m2: 值不对称边下三角
m3: 值不对称边上三角
*/
void My_Matrix_Sys2(int n, int nnzR, int &m1, int &m2, int &m3,
                    std::vector<int>& Row_Number_coo, std::vector<int>& Col_Number_coo, std::vector<VALUE_TYPE>& Val_Number_coo,              
                   std::vector<int>& Row_Number_my_sys_for_fegrass_angle, std::vector<int>& Col_Number_my_sys_for_fegrass_angle, std::vector<VALUE_TYPE>& Val_Number_my_sys_for_fegrass_angle,
                   std::vector<int>& Issystem)
{
    struct timeval t1, t2;
    int i,j,k, p1, p2, p3, tar, tar1, tar2, temp; m1 = m2 = m3 = 0;
    VALUE_TYPE w;
    double time_MY_SYS = 0;
    gettimeofday(&t1, NULL);
    for(i = 0; i < nnzR; ++i)
    {
        if(Issystem[i] == 1 && Row_Number_coo[i] > Col_Number_coo[i]) ++m1;
        else if(Issystem[i] >= 2 && Row_Number_coo[i] > Col_Number_coo[i]) ++m2, ++m3;
    }
    Row_Number_my_sys_for_fegrass_angle.resize(m1+m2+m3);
    Col_Number_my_sys_for_fegrass_angle.resize(m1+m2+m3);
    Val_Number_my_sys_for_fegrass_angle.resize(m1+m2+m3);
    p1 = 0, p2 = m1, p3 = m1 + m2;
    for(int a = 0; a < nnzR; ++a)
    {
        tar = -1;
        if(Row_Number_coo[a] <= Col_Number_coo[a]) continue;
        i = Row_Number_coo[a], j = Col_Number_coo[a], w = Val_Number_coo[a];
        if(Issystem[a] == 1)
        {
            tar = p1++;
            Row_Number_my_sys_for_fegrass_angle[tar] = i;
            Col_Number_my_sys_for_fegrass_angle[tar] = j;
            Val_Number_my_sys_for_fegrass_angle[tar] = w;
        }
        else if(Issystem[a] >= 2)
        {
            tar1 = p2++;
            tar2 = p3++;
            Row_Number_my_sys_for_fegrass_angle[tar1] = i;
            Col_Number_my_sys_for_fegrass_angle[tar1] = j;
            Val_Number_my_sys_for_fegrass_angle[tar1] = w;
            Row_Number_my_sys_for_fegrass_angle[tar2] = j;
            Col_Number_my_sys_for_fegrass_angle[tar2] = i;
            Val_Number_my_sys_for_fegrass_angle[tar2] = Val_Number_coo[Issystem[a]-2];
        }
    }
    gettimeofday(&t2, NULL);
    time_MY_SYS = (t2.tv_sec - t1.tv_sec) * 1000.0 + (t2.tv_usec - t1.tv_usec) / 1000.0;
    // printf("time_MY_SYS_for_fegrass=%lf ms\n", time_MY_SYS);
}

int func1(int n_matrix, int m1, int m2, int m3, int m4, 
    std::vector<int>& Row_Num, std::vector<int>& Col_Num, std::vector<VALUE_TYPE>& Val_Num,
    std::vector<int>& Row_Num_end_fegrass, std::vector<int>& Col_Num_end_fegrass, std::vector<VALUE_TYPE>& Val_Num_end_fegrass,
    double Sratio)
{   
    int n, m, i, i1, i2;
    VALUE_TYPE w;

    n = n_matrix;
    m = Row_Num.size();
	matGraph mat(n, m);
	for(int j = 0; j < m; j++)
	{   
        i1 = Row_Num[j];
        i2 = Col_Num[j];
        w = Val_Num[j];
        mat.addEdge(i1, i2, w);
	}
	mat.sparsify1(Sratio, m1, m2, m3, m4, Row_Num_end_fegrass, Col_Num_end_fegrass, Val_Num_end_fegrass);
	return 0;
}
int func2(int n_matrix, int m1, int m2, int m3, 
    std::vector<int>& Row_Num, std::vector<int>& Col_Num, std::vector<VALUE_TYPE>& Val_Num,
    std::vector<int>& Row_Num_end_fegrass, std::vector<int>& Col_Num_end_fegrass, std::vector<VALUE_TYPE>& Val_Num_end_fegrass,
    double Sratio)
{   
    int n, m, i, i1, i2;
    VALUE_TYPE w;

    n = n_matrix;
    m = Row_Num.size();
	matGraph2 mat(n, m);
	for(int j = 0; j < m; j++)
	{   
        i1 = Row_Num[j];
        i2 = Col_Num[j];
        w = Val_Num[j];
        mat.addEdge(i1, i2, w);
	}
	mat.sparsify1(Sratio, m1, m2, m3, Row_Num_end_fegrass, Col_Num_end_fegrass, Val_Num_end_fegrass);
	return 0;
}

int main(int argc, char **argv)
{   
    char *filename = argv[1];//输入矩阵的路径
    char name[1000];
    int n, m, m1, m2, m3, m4, temp,i;
    double Dtemp;
    double ratio = 0.01;
    std::vector<int> Row_Num_Coo, Col_Num_Coo, Row_Num1, Col_Num1, Row_Num_csp, Col_Num_csp, Row_Num_gpscp, Col_Num_gpscp;
    std::vector<double> Val_Num_Coo, Val_Num1, Val_Num_csp, Val_Num_gpscp;
    // 读取矩阵COO格式
    readMtx(argv[1], n, m, Row_Num_Coo, Col_Num_Coo, Val_Num_Coo);
    getMtxName1(argv[1], name);//zyx
    std::vector<int> Issystem(m, 0);
    getIssystem(n, m, Row_Num_Coo, Col_Num_Coo, Val_Num_Coo, Issystem);

    struct timeval tstart1,tend1;
    gettimeofday(&tstart1, NULL);
    My_Matrix_Sys1(n, m, m1, m2, m3, m4, Row_Num_Coo, Col_Num_Coo, Val_Num_Coo, Row_Num1, Col_Num1, Val_Num1, Issystem);
    // struct timeval tstart1,tend1;
    // gettimeofday(&tstart1, NULL);
    func1(n, m1, m2, m3, m4, Row_Num1, Col_Num1, Val_Num1, Row_Num_csp, Col_Num_csp, Val_Num_csp, ratio);
    gettimeofday(&tend1, NULL);
    double timeUsed_csp = 0;
    timeUsed_csp=1000.0*(tend1.tv_sec-tstart1.tv_sec)+(tend1.tv_usec-tstart1.tv_usec)/1000.0;
    My_Delete1(n, Issystem, Row_Num_Coo, Col_Num_Coo, Val_Num_Coo, Row_Num_csp, Col_Num_csp, Val_Num_csp);
    Row_Num1.clear();
    Col_Num1.clear();
    Val_Num1.clear();

    struct timeval tstart2,tend2;
    gettimeofday(&tstart2, NULL);
    My_Matrix_Sys2(n, m, m1, m2, m3, Row_Num_Coo, Col_Num_Coo, Val_Num_Coo, Row_Num1, Col_Num1, Val_Num1, Issystem);
    // struct timeval tstart2,tend2;
    // gettimeofday(&tstart2, NULL);
    func2(n, m1, m2, m3, Row_Num1, Col_Num1, Val_Num1, Row_Num_gpscp, Col_Num_gpscp, Val_Num_gpscp, ratio);
    gettimeofday(&tend2, NULL);
    double timeUsed_gpscp = 0;  
    timeUsed_gpscp=1000.0*(tend2.tv_sec-tstart2.tv_sec)+(tend2.tv_usec-tstart2.tv_usec)/1000.0;
    My_Delete2(n, Issystem, Row_Num_Coo, Col_Num_Coo, Val_Num_Coo, Row_Num_gpscp, Col_Num_gpscp, Val_Num_gpscp);
    Row_Num1.clear();
    Col_Num1.clear();
    Val_Num1.clear();
    cout << Col_Num_Coo.size() << ' ' << Col_Num_csp.size() << ' ' << Col_Num_gpscp.size() << endl;
    // cout << timeUsed_gpscp << " " << timeUsed_csp << endl;
    saveMtx(saveMtxPath1, n, Row_Num_Coo, Col_Num_Coo, Val_Num_Coo);
    saveMtx(saveMtxPath2, n, Row_Num_csp, Col_Num_csp, Val_Num_csp);
    saveMtx(saveMtxPath3, n, Row_Num_gpscp, Col_Num_gpscp, Val_Num_gpscp);
    saveMtx1(saveMtxPath4, name);//zyx

    string filepath = "/home/zyx/ZYX_iccad24/getData/getdata.csv";
    freopen(filepath.c_str(), "a", stdout);
    char* str = filename;
    printf("%s,gpscp_time = %lf,csp_time = %lf\n",str,timeUsed_gpscp,timeUsed_csp);
    fclose(stdout);

    Row_Num_Coo.clear();
    Col_Num_Coo.clear();
    Val_Num_Coo.clear();
    Row_Num_csp.clear();
    Col_Num_csp.clear();
    Val_Num_csp.clear();
    Row_Num_gpscp.clear();
    Col_Num_gpscp.clear();
    Val_Num_gpscp.clear();
    Issystem.clear();
    return 0;
}
