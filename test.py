import torch
import matplotlib.pyplot as plt

# 定义训练数据
X = torch.tensor([[1.0, 2.0], [3.0, 4.0]],requires_grad=True)
Y = torch.tensor([[2, 4], [6, 8]])


Y_pred=torch.zeros(4)
# 定义损失函数
loss_fn = torch.nn.MSELoss()

# 定义优化器和超参数
optimizer = torch.optim.SGD([X], lr=0.01)
epochs = 100


# 批量梯度下降
for epoch in range(epochs):
    Z=torch.tensor([[2.0, 4.0], [6.0, 8.0]])
    Z[0][0]=X[0][0]+X[1][1]
    loss=torch.pow(torch.det(torch.inverse(X)@torch.inverse(X)-Y),2)+torch.det(Z)
    loss.backward()

