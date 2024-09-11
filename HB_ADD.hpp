#include <vector>

void My_Delete1(int n, std::vector<int>& Issystem,
            std::vector<int> &Row_Number_coo, std::vector<int> &Col_Number_coo, std::vector<double> &Val_Number_coo,
            std::vector<int> &Row_Number_my_sys_fegrass_add_end, std::vector<int>& Col_Number_my_sys_fegrass_add_end, std::vector<double>& Val_Number_my_sys_fegrass_add_end)
{
    int i,j,chosen;
    int nnzR2 = Row_Number_coo.size();
    for(i = 0; i < nnzR2; i++)
    {
        if(Row_Number_coo[i] == Col_Number_coo[i])
        {
            Row_Number_my_sys_fegrass_add_end.push_back(Row_Number_coo[i]);
            Col_Number_my_sys_fegrass_add_end.push_back(Col_Number_coo[i]);
            Val_Number_my_sys_fegrass_add_end.push_back(Val_Number_coo[i]);
        }
    }
}

void My_Delete2(int n, std::vector<int>& Issystem,
            std::vector<int> &Row_Number_coo, std::vector<int> &Col_Number_coo, std::vector<double> &Val_Number_coo,
            std::vector<int> &Row_Number_my_sys_fegrass_add_end, std::vector<int>& Col_Number_my_sys_fegrass_add_end, std::vector<double>& Val_Number_my_sys_fegrass_add_end)
{
    int i,j,chosen;
    int nnzR2 = Row_Number_coo.size();
    for(i = 0; i < nnzR2; i++)
    {
        if(Issystem[i] == 0 || Row_Number_coo[i] == Col_Number_coo[i])
        {
            Row_Number_my_sys_fegrass_add_end.push_back(Row_Number_coo[i]);
            Col_Number_my_sys_fegrass_add_end.push_back(Col_Number_coo[i]);
            Val_Number_my_sys_fegrass_add_end.push_back(Val_Number_coo[i]);
        }
    }
}