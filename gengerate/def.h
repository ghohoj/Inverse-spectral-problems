/*
 * @Description: 一些变量的定义
 * @Author: catgod
 * @Date: 2023-09-05 21:34:45
 * @LastEditTime: 2023-09-06 11:23:41
 * @FilePath: /Inverse spectral problems/gengerate/def.h
 */
#pragma once
#include<fstream>
#include<iostream>
#include<vector>
#include<iomanip>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;
typedef double Real;

//把0到1划分为多少个点
const int precise=5;
//矩阵
typedef Matrix<Real,precise,precise> mat;
//向量
typedef Vector<Real,precise> vec;