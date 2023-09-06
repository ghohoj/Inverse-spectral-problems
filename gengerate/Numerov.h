/*
 * @Description: 实现对 y''+q(x)y=\lambda y y(0)=0 y(1)=0,
    参考论文:Matrix Numerov method for solving Schrödinger’s equation
 * @Author: catgod
 * @Date: 2023-09-02 11:39:30
 * @LastEditTime: 2023-09-06 12:00:36
 * @FilePath: /Inverse spectral problems/gengerate/Numerov.h
 */

#pragma once
#include"def.h"

vec randomvec(){
    vec tmp;
    for(int i=0;i<precise;i++){
        tmp(i)=random();
    }
    return tmp.normalized();
}

/// @brief 求解下三角矩阵的方程组的方程组，Ly=u
/// @param L 下三角矩阵
/// @param u Ly=u中的u
/// @return 返回的y
vec LSolve(const mat&L,const vec& u){
    vec y;
    Real tmp;
    for(int i=0;i<precise;i++){
        tmp=u(i);
        for(int j=0;j<i;j++){
            tmp-=L(i,j)*y(j);
        }
        y(i)=u(i)/L(i,i);
    }
    return y;
}


/// @brief 求解上三角矩阵的方程组，Rx=y
/// @param R 上三角矩阵
/// @param y Rx=y中的y
/// @return 返回的x
vec RSolve(const mat&R,const vec& y){
    vec x;
    Real tmp;
    for(int i=precise-1;i>=0;i--){
        tmp=y(i);
        for(int j=0;j<precise-1-i;j++){
            tmp-=R(i,precise-1-j)*y(j);
        }
        x(i)=y(i)/R(i,i);
    }
    return x;
}

/// @brief 求解Mx=u -> LUx=u -> Ly=u, Ux=y
/// @param L 下三角矩阵
/// @param U 上三角矩阵
/// @param u 求解Mx=u中的u
/// @return 求解的结果x
vec InverseSolve(const mat& L,const mat& U,const VectorXd& u){
    vec y,x;
    y=LSolve(L,u);
    x=RSolve(U,y);
    return x;
}



/// @brief LU分解，即找到M=LU
/// @param M 原来的矩阵
/// @param L 下三角矩阵
/// @param U 上三角矩阵
void LU(const mat& M,mat& L,mat& U){
    Eigen::PartialPivLU<mat> lu(M);
    mat l = mat::Identity();
    L= lu.matrixLU().triangularView<Lower>();
    U = lu.matrixLU().triangularView<Upper>();
}

/// @brief 采取反幂法求解特征值(详细数学推理见readme)
/// @param M 矩阵
/// @param x 需要多少个特征值向量M.rows()
/// @param result 返回的特征值
void InversePowerMethod(const mat& M,int evalnum,vector<Real>& eigval){
    int i=0;
    vector<vec> eigvector;
    vec tmp;
    vec u;
    Real lambda;
    mat L,U;
    LU(M,L,U);
    while (i<evalnum)
    {

        u=randomvec();
        tmp=VectorXd::Zero(M.rows());
        while ((tmp-u).norm()>10e-12)
        {
            tmp=u;
            for(auto v:eigvector){
                u=(u-v.dot(u)*v).normalized();
            }
            if(u(0)<0){
                u=-u;
            }
            u=InverseSolve(L,U,u).normalized();
        }
        lambda=u(0)/InverseSolve(L,U,u)(0);
        eigvector.push_back(u);
        eigval.push_back(lambda);
        i++;
    }   
}


/// @brief see "Matrix Numerov method for solving Schrödinger’s equation"
/// @param q -y''+q(x)y=l y y(0)=0 y(1)=0,q是对q(x)的采样,必须是等间隔的,头尾点必须都采样到
/// @param result 返回值:特征值 
/// @param evalnum 所求特征值的数目
void MatrixNumerov(const vec& q,vector<Real>& eigval,const int& evalnum=10){
    mat Id = mat::Identity();
    mat One1=mat::Zero();
    for(int i=0;i<precise-1;i++){
        One1(i,i+1)=1;
    }
    mat V=mat::Zero();
    for(int i=0;i<precise;i++){
        V(i,i)=q(i);
    }
    auto A=(10*Id+One1+One1.transpose())/12;
    auto B=(-2*Id+One1+One1.transpose())/(precise-1)/(precise-1);
    auto H=B.inverse()*A+V;
    cout<<H<<endl;
    EigenSolver<mat> es(H);
    MatrixXd value= es.eigenvalues().real();
    MatrixXd vector= es.eigenvectors().real();
    cout << value<< endl << endl;
	cout << vector<< endl << endl;
    InversePowerMethod(H,evalnum,eigval);
}




