import numpy as np

from utils.io import *
from utils.scdg import *
import os
import torch
import torch.utils.data as Data
import torch.nn as nn
import torch.nn.functional as F
from sklearn.preprocessing import MinMaxScaler

sptFile = "h5ads/FFPE_Mouse_Kidney.h5spt"
scdata = Load_spt_sc_to_AnnData(sptFile)
stdata = Load_spt_to_AnnData(sptFile, count="SPCS_mat")
scgenes = set(scdata.var_names)
stgenes = set(stdata.var_names)
genes = list(stgenes.intersection(scgenes))
scRef, label = Make_References(sptFile, use_genes=genes)
scRef = np.log(scRef + 1)
count_size = len(genes)
h_dim = 400
z_dim = 20
learning_rate = 1e-3
stdata = stdata[:, genes]
scRef = scRef[genes]
sc.pp.log1p(stdata)
stdataset = stdata.X.todense()

# DCV model
class DCV(nn.Module):
    def __init__(self, count_size = 3000, h_dim=400, z_dim=20, o_dim=350):
        super(DCV, self).__init__()
        self.fc1 = nn.Linear(count_size, h_dim)
        self.fc2 = nn.Linear(h_dim, z_dim)  # 均值 向量
        self.fc3 = nn.Linear(h_dim, z_dim)  # 保准方差 向量
        self.fc4 = nn.Linear(z_dim, h_dim)
        self.fc5 = nn.Linear(h_dim, o_dim)

    # 编码过程
    def encode(self, x):
        h = F.relu(self.fc1(x))
        return self.fc2(h), self.fc3(h)

    # 随机生成隐含向量
    def reparameterize(self, mu, log_var):
        std = torch.exp(log_var / 2)
        eps = torch.randn_like(std)
        return mu + eps * std

    # 解码过程
    def decode(self, z):
        h = F.relu(self.fc4(z))
        return torch.softmax(self.fc5(h), dim=1)

    # 整个前向传播过程：编码-》解码
    def forward(self, x):
        mu, log_var = self.encode(x)
        z = self.reparameterize(mu, log_var)
        x_reconst = self.decode(z)
        return x_reconst, mu, log_var

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
num_epoches = 1000
epoch_thresh = 1e-5
uni_types = np.unique(stdata.obs['BayesSpace'])
for type in uni_types:
    model = DCV(count_size=count_size, h_dim=h_dim, z_dim=z_dim, o_dim=len(scRef)).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)
    zidx = np.where(stdata.obs['BayesSpace'] == type)[0]
    dataset = stdataset[zidx, :]
    dataset = torch.from_numpy(np.array(dataset))
    ref = torch.from_numpy(np.array(scRef, dtype="float32"))
    batch_size = len(dataset)
    data_loader = Data.DataLoader(dataset=dataset,
                                  batch_size=batch_size,
                                  shuffle=True)
    for epoch in range(num_epoches):
        for i, x in enumerate(data_loader):
            # 获取样本，并前向传播
            x = x.to(device).view(-1, count_size)
            x_reconst, mu, log_var = model(x)
            reconst_loss = F.mse_loss(torch.matmul(x_reconst, ref), x, reduction='mean')
            kl_div = - 0.5 * torch.sum(1 + log_var - mu.pow(2) - log_var.exp())
            # 反向传播和优化
            loss = reconst_loss + kl_div
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

            if epoch % 100 == 0 and epoch > 0:
                print("Epoch[{}/{}] Reconst Loss: {}, KL Div: {}"
                      .format(epoch, num_epoches, reconst_loss.item(), kl_div.item()))
    print("Epoch has reach an end.")
    with torch.no_grad():
        # 重构的图像
        out, _, _ = model(x)
        npout = out.numpy()
        probs = []
        for i in range(0, npout.shape[1], 25):
            prob = np.sum(npout[:, i:i+25], axis=1)
            probs.append(prob.T)
        probs = np.array(probs).T

        reconst = torch.matmul(out, ref)

