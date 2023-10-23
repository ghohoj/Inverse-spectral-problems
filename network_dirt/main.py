'''
Description: 展示直接做泰勒展开的结果
Author: catgod
Date: 2023-09-08 15:18:22
LastEditTime: 2023-10-21 23:17:47
FilePath: /Inverse spectral problems/network_dirt/main.py
'''

import numpy as np
import torch
import  torch.nn as nn
import  torch.optim as optim
import torch.utils.data as Data
import json


train_data=[]
target_data=[]
for j in range(10):
    with open('./data/test'+str(j)+'.json') as f:
        loaddata = json.load(f)
        for i in loaddata:
            train_data.append(i["lambda"])
            target_data.append(i["cof"])
train_data = torch.tensor(train_data,dtype=torch.float)
target_data = torch.tensor(target_data,dtype=torch.float)
torch_dataset = Data.TensorDataset(train_data, target_data)


with open('./data/test.json') as f:
    loaddata = json.load(f)
test_data=[]
test_target_data=[]
for i in loaddata:
    test_data.append(i["lambda"])
    test_target_data.append(i["cof"])
test_data = torch.tensor(test_data,dtype=torch.float)
test_target_data = torch.tensor(test_target_data,dtype=torch.float)
test_dataset = Data.TensorDataset(test_data, test_target_data)


data_num=np.shape(train_data)[0]
feature_num=np.shape(train_data)[1]
target_num=np.shape(target_data)[1]
normal=torch.tensor([((i+1)*np.pi)**2 for i in range(feature_num)],dtype=torch.float)
normalsize=torch.tensor([(1/(i+1))**2 for i in range(target_num)],dtype=torch.float)
learning_rate=0.005
epochs=100
rangenum=10
BATCH_SIZE=16


train_loader = Data.DataLoader(
    dataset=torch_dataset,      # torch TensorDataset format
    batch_size=BATCH_SIZE,      # mini batch size
    # shuffle=True,               # 要不要打乱数据 (打乱比较好)
    num_workers=2,              # 多线程来读数据
)


class MLP(nn.Module):#自定义类 继承nn.Module
    def __init__(self):#初始化函数
        super(MLP, self).__init__()#继承父类初始化函数
        self.model1 = nn.Sequential(
            nn.Linear(feature_num, feature_num),
            nn.ELU(),
            nn.Linear(feature_num, feature_num),
            nn.ELU(),
            nn.Linear(feature_num, target_num),
        )
        self.model2 = nn.Sequential(
            nn.Linear(feature_num, feature_num),
            nn.ELU(),
            nn.Linear(feature_num, feature_num),
            nn.ELU(),
            nn.Linear(feature_num, target_num),
        )
    def forward(self, x):
        tmp=x-normal
        y = self.model1(tmp)
        y = (y+self.model2(tmp))
        return y

""" def losscal(coff:torch.Tensor,target:torch.Tensor):
    result=torch.norm(coff-target)
    return result """

criterion = torch.nn.MSELoss(size_average=False)

net = MLP()#创建实例

print(net)
optimizer = optim.Adam(net.parameters(), lr=learning_rate)
for epoch in range(epochs):
    for batch_idx, (data, target_data) in enumerate(train_loader):
        coefficient = net(data)#得到x经过模型后的输出
        train_loss =criterion(coefficient,target_data)
        optimizer.zero_grad()
        train_loss.backward()
        optimizer.step()
    print('epoch',epoch,float(train_loss.data))


print(criterion(net(test_data),test_target_data))