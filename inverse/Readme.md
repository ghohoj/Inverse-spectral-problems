# 问题

已知特征值l，
求解
$$
-y''+q(x)y=l y\\ y(0)=0 , y(1)=0
$$
中的q

# 代码思路

无，采用Numerical algorithms for inverse Sturm-Liouville problems复现即可。

# 测试

## 失败测试

1. pytorch禁止自替代运算（内存上重复的算符，例如a=inv(a)），例如 x=f(q)之后又写 x=g(q)，不应该出现+=等自替代符号，以下代码也是违法的
    ```
    self.TrainPower=torch.zeros((s,K1,K1))
    self.TrainPower[i]=self.TrainPower[i]+torch.mm(self.TrainPower[1],self.TrainPower[i-1])
    ```

    原因是前后出现自替代,改为

    ```
    self.TrainPower=[torch.zeros((K1,K1)) for i in range(s+1)]
    self.TrainPower[i]=self.TrainPower[i]+torch.mm(self.TrainPower[1],self.TrainPower[i-1])
    ```

## 成功测试

假定特征值为$(n\pi)^2$，众所周知q=0，测试运行test.py/test1()，成功
