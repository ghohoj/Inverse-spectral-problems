'''
Description: 论文Numerical algorithms for inverse Sturm-Liouville problems的复现。
Author: catgod
Date: 2023-09-05 21:44:27
LastEditTime: 2023-09-09 22:06:28
FilePath: /Inverse spectral problems/inverse/InverseSLProblem.py
'''
from typing import List, Dict
import numpy as np
import torch
from torchmin import minimize
import json
import matplotlib.pyplot as plt 


torch.set_printoptions(precision=8)
N=1

class PPoly:
    def __init__(self,s:int,K1:int) -> None:
        """论文Numerical algorithms for inverse Sturm-Liouville problems(2.5)(2.6)公式,主要功能为实现(2.7)两边的计算
        
        Args:
            s (int): P的最高阶数
            K1 (int): 特征值的数目
        """
        self.s=s
        self.K1=K1
        # P多项式的系数，P[i][j]为第i个多项式子的x^j系数
        self.PPar=torch.zeros((s,s+1))
        # TargetPower[i][j] 第i个特征值的j次方
        self.TargetPower=torch.zeros((K1,s+1))
        # TrainPower[i] 矩阵的i次方
        self.TrainPower=[torch.zeros((K1,K1)) for i in range(s+1)]
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
        for i in range(self.s+1):
            if i==0:
                self.TrainPower[i]=self.TrainPower[i]+torch.eye(np.shape(self.TrainPower[-1])[0])
            elif i==1:
                self.TrainPower[i]=self.TrainPower[i]+torch.inverse(Mat)
            else:
                self.TrainPower[i]=self.TrainPower[i]+self.TrainPower[1]@self.TrainPower[i-1]
     
    def polyEval(self)->torch.Tensor:
        """计算r_true,定义见论文Numerical algorithms for inverse Sturm-Liouville problems Algorithm 1第四行

        Returns:
            torch.Tensor: 计算r_true
        """
        result=torch.zeros((self.s))
        for j in range(self.s):
            for x in self.TargetPower:
                for k in range(self.s+1):
                    result[j]=result[j]+x[k]*self.PPar[j][k]
        return result
    def polyMAT(self)->torch.Tensor:
        """计算r,定义见论文Algorithm 1第17行

        Returns:
            torch.Tensor: 计算r
        """
        result=torch.zeros((self.s))
        tmpresult=torch.zeros((self.s,
                            self.TrainPower[0].size()[0],
                            self.TrainPower[0].size()[1]))
        for j in range(self.s):
            for k in range(len(self.TrainPower)):
                tmpresult[j]=tmpresult[j]+self.TrainPower[k]*self.PPar[j][k]
            result[j]=result[j]+torch.trace(tmpresult[j])
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


def read(address:str="data.json")->Dict[str,List[float]]:
    """读取json文件

    Args:
        address (str): 文件地址

    Returns:
        Dict[str,List[float]]: {q:(0,1)中q的值,lambda:特征值}
    """
    with open('./data/'+address, 'r') as fcc_file:
        fcc_data = json.load(fcc_file)
    return fcc_data


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

    D=torch.zeros((P.K1,P.K1))
    for i in range(P.K1):
        D[i][i]=D[i][i]+((i+1)*np.pi)**2
    M=torch.zeros((P.K1,P.K1))
    for i in range(P.K1):
        for j in range(P.K1):
            for k in range(N):
                M[i][j]=M[i][j]+q[k]*MIntegral(i+1,j+1,k+1)
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
    result:List[float]=[];
    q_result=q_result.detach().numpy();
    for i in range(precise):
        tmp=0
        for j in range(np.size(q_result)):
            tmp+=q_result[j]*np.cos(2*j*np.pi*(1.0*i/precise))
        result.append(tmp)
    return result

def train(x:torch.Tensor,q_num:int=N):
    """训练函数

    Args:
        x (torch.Tensor): 特征值
        q_num (int):展开的项数
    Returns:
        _type_: q(x)=sum_{k=1} q_k cos(2(k-1)pi x)中的q_k
    """
    q=torch.ones((q_num))*(1)
    P=PPoly(len(x),len(x))
    target=gettarget(x,P)
    #quasi-newton,调用求解器
    # q_result = minimize(lambda q:forward(q,target,P), q, method='l-bfgs',
    #                     options=dict(line_search='strong-wolfe'), tol=1e-6,
    #                     max_iter=50,disp=2)
    q.requires_grad = True
    optimizer = torch.optim.SGD([q], lr=1)
    for epoch in range(200):
        P=PPoly(len(x),len(x))
        loss=forward(q,target,P)
        loss.backward()
        optimizer.step()
        optimizer.zero_grad()
        print("第"+str(epoch)+"轮迭代：")
        print(q)
        print(loss)
        #变梯度下降
        if(epoch==30):
            for ps in optimizer.param_groups:
                ps["lr"]=0.5
        if(epoch==60):
            for ps in optimizer.param_groups:
                ps["lr"]=0.1
        if(epoch==90):
            for ps in optimizer.param_groups:
                ps["lr"]=0.05
        if(epoch==120):
            for ps in optimizer.param_groups:
                ps["lr"]=0.01
        if(epoch==150):
            for ps in optimizer.param_groups:
                ps["lr"]=0.005
    return q




def main():
    """主函数
    """
    data=read()
    q_result=train(torch.tensor(data["lambda"]))
    plt.plot(qktoqx(q_result,3000))
    plt.plot(np.array([0]*3000))
    plt.ylim([-1,5])
    plt.show()
    plt.pause(0)
    
main()