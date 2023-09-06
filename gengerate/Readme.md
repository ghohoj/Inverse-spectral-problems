采取反幂法求解更多特征值

在第一轮中

```
u[i+1] =M^-1 u[i]
```

多次迭代得到第一个特征向量记录为v，特征值为i

在第i轮中，我们需要反复排除前面已经出现的特征，我们假设已经加入的v向量都是正交的

```
u[i]=u[i]-sumforv(v dot u[i])v
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
        if(result is empty){
        	vectorxd u=(1,0,0,...)
        }
        else{
        	vectorxd u=选取与已经加入result的向量垂直的任意向量
        }
        while(norm(tmp-u)>pre){
        	tmp=u;
        	for(v:result){
        		u=norm(u-(v dot u)v)
        	}
        	u=M^-1 u;
        }
        result.push(u的特征值);
        i++;
    }
    
}
```

这个代码就可以从小到大求解特征值与对应的特征向量了。

```
u=M^-1 u
```

这里一方面可以调用库，或在采取LU分解，或在在论文中有更为优化的算法，都可以。

