'''
Description: 论文Numerical algorithms for inverse Sturm-Liouville problems的复现,注意要求关于0.5的轴对称！
Author: catgod
Date: 2023-09-05 21:44:27
LastEditTime: 2023-09-06 10:05:21
FilePath: /Inverse spectral problems/inverse/InverseSLProblem.py
'''
from typing import List, Dict, Tuple
import numpy as np
import torch
from torchmin import minimize
import json


N=10

class PPoly:
    def __init__(self,s:int,K1:int) -> None:
        """论文Numerical algorithms for inverse Sturm-Liouville problems(2.5)(2.6)公式,主要功能为实现(2.7)两边的计算
        
        Args:
            s (int): P的最高阶数
            K1 (int): 特征值的数目
        """
        self.s=s
        # P多项式的系数，P[i][j]为第i个多项式子的x^j系数
        self.PPar=torch.zeros((s,s+1))
        # TargetPower[i][j] 第i个特征值的j次方
        self.TargetPower=torch.zeros((K1,s+1))
        # TrainPower[i] 矩阵的i次方
        self.TrainPower=torch.zeros((s+1,K1,K1))
        self._PPoly()

    def _PPoly(self)->None:
        """计算P多项式从0到s的所有阶的系数,存储在self.PPar[i][j]
        """
        self.PPar[0][1]=1;
        for i in range(1,self.s):
            tmp=(-1)**(i+1)
            for j in range(i+2):
                self.PPar[i][j]=self.PPar[i-1][j-1]+tmp*self.PPar[i-1][j]


    def getTargetPower(self,myeval:List[float])->None:
        """计算特征值的幂(0,-1,-2,...),记录在TargetPower

        Args:
            myeval (List[float]): 特征值(从小到大)
        """
        for i in range(len(myeval)):
            for j in range(self.s):
                if j==0:
                    self.TargetPower[i][j]=1
                elif j==1:
                    self.TargetPower[i][j]=1/myeval[i]
                else:
                    self.TargetPower[i][j]=self.TargetPower[i][j-1]/myeval[i]

    def getTrainPower(self,Mat:torch.Tensor)->None:
        """计算矩阵的幂(0,-1,-2,...),记录在TrainPower

        Args:
            Mat (torch.Tensor): 矩阵R(q)
        """
        for i in range(self.s):
            if i==0:
                self.TrainPower[i]=torch.eye(np.shape(self.TrainPower[-1])[0])
            elif i==1:
                self.TrainPower[i]=torch.linalg.inv(Mat)
            else:
                self.TrainPower[i]=self.TrainPower[1]@self.TrainPower[i-1]
     
    def polyEval(self)->torch.Tensor:
        """计算r_true,定义见论文Numerical algorithms for inverse Sturm-Liouville problems Algorithm 1第四行

        Returns:
            torch.Tensor: 计算r_true
        """
        result=torch.zeros((self.s))
        for j in range(self.s):
            tmpresult=0
            for x in self.TargetPower:
                for k in range(self.s+1):
                    tmpresult+=x[k]*self.PPar[j][k]
            result[j]=tmpresult
        return result
    def polyMAT(self)->torch.Tensor:
        """计算r,定义见论文Algorithm 1第17行

        Returns:
            torch.Tensor: 计算r
        """
        result=torch.zeros((self.s))
        for j in range(self.s):
            tmpresult=torch.zeros_like(self.TrainPower[0])
            for k in range(len(self.TrainPower)):
                tmpresult+=self.TrainPower[k]*self.PPar[j][k]
            result[j]=torch.trace(tmpresult)
        return result


def MIntegral(i:int,j:int,k:int)->float:
    """论文Numerical algorithms for inverse Sturm-Liouville problems(3.11)公式

    Args:
        i (int): i
        j (int): j
        k (int): k

    Returns:
        float: M(e_K)_{ij}
    """
    result=0
    if i+j+2*k-2==0:
        result=-0.5
    elif i+j-2*k+2==0:
        result=-0.5
    elif i-j+2*k-2==0 and k!=1:
        result=0.5
    elif i-j-2*k+2==0 and k!=1:
        result=0.5
    elif i==j and k==1:
        result=1
    else:
        result=0
    return result


def read(address:str)->Dict[str,List[float]]:
    """读取json文件

    Args:
        address (str): 文件地址

    Returns:
        Dict[str,List[float]]: {q:(0,1)中q的值,lambda:特征值}
    """
    return 



def forward(q:torch.Tensor,target:torch.Tensor,P:PPoly)->torch.float:
    """计算优化的目标

    Args:
        q (torch.Tensor): Algorithhm 中的q,q(x)=sum_{k=1} q_k cos(2(k-1)pi x)中的q_k
        target(torch.Tensor): Algorithhm 中的r_true
        P (PPoly):多项式
    Returns:
        torch.float: Algorithhm 中的F(q)=|r-r_{true}|,最小化的目标
    """
    #计算R矩阵(R=D+M,求逆在PPoly进行)
    D=torch.zeros((N,N))
    for i in range(N):
        D[i][i]=((i+1)*np.pi)**2
    M=torch.zeros((N,N))
    for i in range(N):
        for j in range(N):
            for k in range(N):
                M[i][j]+=q[k]*MIntegral(i,j,k)
    R=D+M
    # 在P多项式中计算trace的结果
    P.getTrainPower(R)
    tmp=P.polyMAT()
    return (tmp-target).norm()


def gettarget(x:torch.Tensor,P:PPoly)->torch.float:
    """调用对象P计算r_true

    Args:
        x (torch.Tensor): 特征值
        P (PPoly): PPoly对象

    Returns:
        torch.float: r_true
    """
    P.getTargetPower(x)
    target=P.polyEval()
    return target


def qktoqx(q_result:torch.Tensor,precise:int)->List[float]:
    """计算q(x)=sum_{k=1} q_k cos(2(k-1)pi x),即q在(0,1)分precise份

    Args:
        q_result (torch.Tensor): q_k
        precise(int): 精度
    Returns:
        List[float]: q(x)
    """
    return 

def train(x:torch.Tensor):
    """训练函数

    Args:
        x (torch.Tensor): 特征值

    Returns:
        _type_: 拟合的结果
    """
    q=torch.zeros((9),requires_grad=True)
    P=PPoly(10,2)

    target=gettarget(x,P)


    #quasi-newton,调用求解器
    q_result = minimize(lambda q:forward(q,target,P), q, method='bfgs')

    result=qktoqx(q_result,1000)
    
    return result
