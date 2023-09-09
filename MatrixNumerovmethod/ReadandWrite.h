/*
 * @Description: 写入文件
 * @Author: catgod
 * @Date: 2023-09-06 09:56:59
 * @LastEditTime: 2023-09-07 20:41:33
 * @FilePath: /Inverse spectral problems/MatrixNumerovmethod/ReadandWrite.h
 */

#pragma once
#include"def.h"

/// @brief 写入json文件
/// @param q 势能函数
/// @param lam 特征值
void write(const vec& q,vector<Real>& lam,string address="data.json")
{
    address="./data/"+address;
	fstream f;
	f.open(address,ios::out);
	f<<"{";
    f<<" \"q\":[ ";
    for(int i=0;i<precise-1;i++){
        f<<setprecision(10)<<q.coeff(i)<<",";
    }
    f<<setprecision(10)<<q.coeff(precise-1);
    f<<"],";
    f<<" \"lambda\":[ ";
    sort(lam.begin(),lam.end());
    for(int i=0;i<lam.size()-1;i++){
        f<<setprecision(10)<<lam[i]<<",";
    }
    f<<setprecision(10)<<lam[lam.size()-1];
    f<<"]";
    f<<"}"<<endl;
	f.close();
}
