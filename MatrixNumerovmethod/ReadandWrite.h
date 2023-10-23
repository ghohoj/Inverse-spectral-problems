/*
 * @Description: 写入文件
 * @Author: catgod
 * @Date: 2023-09-06 09:56:59
 * @LastEditTime: 2023-09-15 16:28:27
 * @FilePath: /Inverse spectral problems/MatrixNumerovmethod/ReadandWrite.h
 */

#pragma once
#include"def.h"

/// @brief 写入数据
/// @param f 文件流
/// @param q 势能函数
/// @param lam 特征值
void getin(fstream& f,const vector<Real>& q,vector<Real>& lam){
    f<<"{";
    f<<" \"q\":[ ";
    for(int i=0;i<precise-1;i++){
        f<<setprecision(10)<<q[i]<<",";
    }
    f<<setprecision(10)<<q[precise-1];
    f<<"],";
    f<<" \"lambda\":[ ";
    sort(lam.begin(),lam.end());
    for(int i=0;i<lam.size()-1;i++){
        f<<setprecision(10)<<lam[i]<<",";
    }
    f<<setprecision(10)<<lam[lam.size()-1];
    f<<"]";
    f<<"}"<<endl;
}

/// @brief 写入json文件
/// @param q 势能函数
/// @param lam 特征值
void write(const vector<Real>& q,vector<Real>& lam,string address="data.json")
{
    address="./data/"+address;
	fstream f;
	f.open(address,ios::out);
	getin(f,q,lam);
	f.close();
}



void write(const vector<vector<Real>>& qs,vector<vector<Real>>& lams,string address="data.json")
{
    address="./data/"+address;
	fstream f;
	f.open(address,ios::out);
    f<<"["<<"\n";
    for(int i=0;i<qs.size();i++){
        getin(f,qs[i],lams[i]);
        if(i!=qs.size()-1){
            f<<","<<"\n";
        }
    }
    f<<"]"<<"\n";
	f.close();
}
