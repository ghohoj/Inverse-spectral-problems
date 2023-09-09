# 问题声明

求解
$$
-y''+q(x)y=l y\\ y(0)=0 , y(1)=0
$$
中的特征值l

# 论文更正声明

基本上采用Numerov method for solving Schrödinger’s equation的方法。

Matrix Numerov method for solving Schrodinger’s equation文章中$B^{-1}$计算量非常夸张，文章中，对一个300*300的矩阵的逆，我无法理解为什么可以他在0.1s内算出来。

一个比较简单的更正是：直接忽略$B^{-1}$。

# 采取反幂法求解更多特征值

在第一轮中，反复求解

```
M u[i+1] = u[i]
```

多次迭代得到第一个特征向量记录为v，特征值为i

在第i轮中，我们需要反复排除前面已经出现的特征，我们假设已经加入的v向量都是正交的

```
u[i]=u[i]-sumforv((v dot u[i]) v)
u[i+1] =M^-1 u[i]
```

以此来排除原先分量的影响

最后得到的向量自然也和已经记录的向量正交

由此我们的伪代码为

```
void InversePowerMethod(const MatrixXd& M,int x,vector<Real>& result){
    int i=0;
    vectorxd tmp=(0,0,0,...);//记录上一轮的u
    while (i<x)
    {
        vectorxd u=随机初始化
        while(norm(tmp-u)>pre){
        	tmp=u;
        	for(v:result){
        		u=norm(u-(v dot u)v)
        	}
        	求解M y= u;
            u=y
        }
        result.push(u的特征值);
        i++;
    }
    
}
```

这个代码就可以从小到大求解特征值与对应的特征向量了。

# 测试

## 失败测试

1. 全程必须采用稀疏矩阵！否则eigen会警告内存不足。
2. 不要尝试以下思路：在第一轮计算之后，运行M+=一个大值*u^T u，以此来扩大原来最小的特征值，进而求出第二小的特征值，这样的问题是，会把原先的稀疏矩阵破坏为一个满矩阵，然后就会算不出来了。
3. 不要尝试求取太多的特征值，原因是收敛的速度与$ \lambda_1/\lambda_2 $的幂有关，要是$\lambda$太大可能会造成$ \lambda_1/\lambda_2 \to 1$，然后就不收敛了。

## 成功测试

运行test1()中，我们让q=0，显然$\lambda=(n\pi)^2$
程序计算结果为
```
迭代:8
迭代:27
迭代:19
迭代:43
9.85645
39.4258 
88.708
157.703
```

基本吻合结果。
