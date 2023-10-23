'''
Description: 请加入备注来说明这个文件的意义
Author: catgod
Date: 2023-10-20 23:30:39
LastEditTime: 2023-10-20 23:30:40
FilePath: /Inverse spectral problems/network_dirt/test.py
'''
import torch

#loss  必须是一个标量才可以用backward。

#1.准备数据集
x_data = torch.Tensor([[1.0],[2.0],[3.0]])
y_data = torch.Tensor([[2.0],[4.0],[6.0]])

#2.使用Class设计模型
class LinearModel(torch.nn.Module):
    #构造函数，初始化
    def __init__(self):
        super(LinearModel,self).__init__()
        self.liner = torch.nn.Linear(1,1)  #(1,1)表示输入和输出的维度
        
    #前向传播函数
    #forward（）相当于对父类_init_（）进行重载   
    def forward(self,x):
        y_pred = self.liner(x)
        return y_pred
    
model = LinearModel()  #创建类LinearModel的实例

#3.构建损失函数和优化器的选择
criterion = torch.nn.MSELoss(size_average=False)
optimizer = torch.optim.SGD(model.parameters(),lr=0.01)

#4.进行训练迭代
for epoch in range(1000):
    y_pred = model(x_data)
    loss = criterion(y_pred,y_data)
    print(epoch,loss.item())  #.item()获取数值大小
    
    #由.backward（）计算的grad将被累积。 因此，在反向传播之前，记住将梯度设置为0。
    optimizer.zero_grad()
    loss.backward()    
    optimizer.step()  #进行更新update
    
#输出权重w和偏置b
print('w=',model.liner.weight.item())
print('b=',model.liner.bias.item())

#测试模型
x_test = torch.Tensor([4.0])
y_test = model(x_test)
print('y_pred = ',y_test.data)
