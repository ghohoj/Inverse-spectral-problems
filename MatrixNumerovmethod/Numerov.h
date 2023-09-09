/*
 * @Description: 实现对 -y''+q(x)y=\lambda y y(0)=0 y(1)=0,
    参考论文:Matrix Numerov method for solving Schrödinger’s equation
 * @Author: catgod
 * @Date: 2023-09-02 11:39:30
 * @LastEditTime: 2023-09-09 17:08:03
 * @FilePath: /Inverse spectral problems/MatrixNumerovmethod/Numerov.h
 */

#pragma once
#include"def.h"


double sum(vector<double>& xs){
    double res=0;
    for(auto x:xs){
        res+=x;
    }
    return res;
} 


/// @brief 采取反幂法求解特征值(详细数学推理见readme)
/// @param M 矩阵
/// @param x 需要多少个特征值向量M.rows()
/// @param result 返回的特征值
void InversePowerMethod(const mat& M,int evalnum,vector<Real>& eigval){
    int i=0;
    vector<vec> eigvector;
    vec u(precise);
    Real lambda;
    mat tmpmat(precise,precise);
    mat tmpmat2(precise,precise);
    while (i<evalnum)
    {
        u=randomvec();
        vec tmp(precise);
        SparseLU<mat> solver;
        mat Id(precise,precise);
        double tmpint;
        Id.setIdentity();
        // if(eigval.empty()){
        //     tmpmat=M;
        //     tmpmat2=tmpmat;
        // }
        // else{
        //     tmpmat=tmpmat+
        //     (10000*eigval[i-1])*eigvector[i-1]*eigvector[i-1].transpose();
        //     tmpmat2=tmpmat-eigval[i-1]*Id;
        // }

        tmpmat2=M;
        solver.compute(tmpmat2);
        
        int cyctimes=0;
        while ((tmp-u).norm()>1.0/precise/1000)
        {
            tmp=u;
            for(auto v:eigvector){
                u=(u-v.dot(u)*v);
                u=u/u.norm();
            }
            if(u.coeff(0)<0){
                u=-u;
            }
            u=solver.solve(u);
            u=u/u.norm();
            cyctimes++;
        }
        cout<<"迭代:"<<cyctimes<<endl;
        tmp=solver.solve(u);
        lambda=u.norm()/tmp.norm();
        eigvector.push_back(u);
        
        // if(eigval.empty()){
        //     tmpint=lambda;
        // }
        // else{
        //     tmpint=lambda+eigval[i-1];
        // }
        // eigval.push_back(tmpint);

        eigval.push_back(lambda);
        i++;
    }   
}

mat Inverse_BA(const mat& B,const mat& A){
    SparseLU<SparseMatrix<double> > solver;
    solver.compute(B);
    auto inv = solver.solve(A);
    return inv;
}


/// @brief see "Matrix Numerov method for solving Schrödinger’s equation"
/// @param q -y''+q(x)y=l y y(0)=0 y(1)=0,q是对q(x)的采样,必须是等间隔的,头尾点必须都采样到
/// @param result 返回值:特征值 
/// @param evalnum 所求特征值的数目
void MatrixNumerov(const vector<Real>& q,vector<Real>& eigval,const int& evalnum=10){
    mat Id(precise,precise);
    Id.setIdentity();
    mat One1(precise,precise);
    for(int i=0;i<precise-1;i++){
        One1.insert(i,i+1)=1;
        One1.insert(i+1,i)=1;
    }
    mat V(precise,precise);
    for(int i=0;i<precise;i++){
        V.insert(i,i)=q[i];
    }
    mat A(precise,precise);
    mat B(precise,precise);
    // A=(10*Id+One1)/12;
    B=(2*Id-One1)*(precise-1)*(precise-1);
    // cout<<B;
    mat H(precise,precise);
    H=B+V;
    InversePowerMethod(H,evalnum,eigval);
}
