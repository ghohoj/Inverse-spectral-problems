/*
 * @Description: 一些变量的定义
 * @Author: catgod
 * @Date: 2023-09-05 21:34:45
 * @LastEditTime: 2023-09-07 15:40:57
 * @FilePath: /Inverse spectral problems/MatrixNumerovmethod/def.h
 */
#pragma once
#include<fstream>
#include<iostream>
#include<vector>
#include<iomanip>
#include<math.h>
#include<random>
#include <Eigen/Dense>
#include<Eigen/SparseLU>
using namespace std;
using namespace Eigen;
typedef double Real;

//把0到1划分为多少个点
const int precise=3000;
//矩阵
typedef SparseMatrix<Real> mat;
//向量
typedef SparseVector<Real> vec;


/// @brief 返回0-1随机数，固定随机种子
/// @return 随机数
double myrand(){
    int x=rand()%100000;
    return x*1.0/100000;
}

/// @brief 返回一个随机向量
/// @return 一个随机向量
vec randomvec(){
    vec tmp(precise);
    for(int i=0;i<precise;i++){
        tmp.insert(i)=myrand();
    }
    return tmp/tmp.norm();
}