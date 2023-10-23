/*
 * @Description: 生成随机的势能
 * @Author: catgod
 * @Date: 2023-09-05 21:31:55
 * @LastEditTime: 2023-09-15 22:36:16
 * @FilePath: /Inverse spectral problems/MatrixNumerovmethod/rV.h
 */

#pragma once
#include"def.h"
#include<random>

/// @brief 多项式计算 a0+a1 x+a2 x^2+a3 x^3+...
/// @param cof (a0,a1,a2,...)
/// @param x x
/// @return 计算结果
Real PolyCal(vector<Real>& cof,Real x){
    auto s=cof.size();
    Real m=1;
    Real result=0;
    for(int i=0;i<s;i++){
        result+=cof[i]*m;
        m*=x;
    }
    return result;
}


/// @brief 随机返回一个势能，采取的是多项式方式(不需要保证函数一直单调或保证函数为正数)
/// @param n 多项式最高n-1阶
/// @return 势能函数，(0,1)之间分为precise份
vector<Real> ployV(int n=4){
    vector<Real> result(precise);
    vector<Real> cof(n);
    for(int i=0;i<n;i++){
        cof[i]=(myrand()-0.5)*20;
    }
    for(int i=0;i<precise;i++){
        result[i]=PolyCal(cof,i*1.0/precise);
    }
    return result;
}

vector<Real> ployV(vector<Real>& cof){
    vector<Real> result(precise);
    for(int i=0;i<precise;i++){
        result[i]=PolyCal(cof,i*1.0/precise);
    }
    return result;
}

/// @brief 返回全为0，从数学上最后生成的特征值应该与n^2成正比
/// @return 势能函数，(0,1)之间分为precise份
vector<Real> testV(){
    vector<Real> tmp(precise);
    for(long int i=0;i<precise;i++){
        tmp[i]=0;
    }
    return tmp;
}


Real CosCal(vector<Real>& cof,Real x){
    auto s=cof.size();
    Real result=0;
    for(int i=0;i<s;i++){
        result+=cof[i]*cos(2*M_PI*i*x);
    }
    return result;
}



vector<Real> CosV(int n=4){
    vector<Real> tmp;
    vector<Real> cof(n);
    for(int i=0;i<n;i++){
        cof[i]=random()-0.5;
    }
    for(int i=0;i<precise;i++){
        tmp[i]=tmp[i]+CosCal(cof,i*1.0/1000);
    }
    return tmp;
}

vector<Real> CosV(vector<Real> cof){
    vector<Real>tmp(precise);
    for(int i=0;i<precise;i++){
        tmp[i]=tmp[i]+CosCal(cof,i*1.0/1000);
    }
    return tmp;
}