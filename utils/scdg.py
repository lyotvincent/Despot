from utils.io import *
import os
import torch
import torch.utils.data as Data
import torch.nn as nn
import torch.nn.functional as F
from sklearn.preprocessing import MinMaxScaler
# Single Cell Referenecs Down Sampling Through Graph Convolution Network


# VAE model
class VAE(nn.Module):
    def __init__(self, count_size = 3000, h_dim=400, z_dim=20):
        super(VAE, self).__init__()
        self.fc1 = nn.Linear(count_size, h_dim)
        self.fc2 = nn.Linear(h_dim, z_dim)  # 均值 向量
        self.fc3 = nn.Linear(h_dim, z_dim)  # 保准方差 向量
        self.fc4 = nn.Linear(z_dim, h_dim)
        self.fc5 = nn.Linear(h_dim, count_size)

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
        return torch.sigmoid(self.fc5(h))

    # 整个前向传播过程：编码-》解码
    def forward(self, x):
        mu, log_var = self.encode(x)
        z = self.reparameterize(mu, log_var)
        x_reconst = self.decode(z)
        return x_reconst, mu, log_var


def Make_References(sptFile, count_size=3000,
                    h_dim=400, z_dim=20,
                    num_epochs=1000, learning_rate=1e-3,
                    epoch_thresh=1e-5, ref_size=25, use_genes=None):
    # 设备配置
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    scdata = Load_spt_sc_to_AnnData(sptFile)
    if use_genes is not None:
        scdata = scdata[:, use_genes]
        count_size = len(use_genes)
    else:
        scdata = scdata[:, scdata.uns['HVGs']]
    sc.pp.log1p(scdata)
    scdataset = scdata.X.todense()
    uni_types = np.unique(scdata.obs['annotation'])
    scaler = MinMaxScaler()
    scdataset = scaler.fit_transform(scdataset)  # 注意，这里的values是array

    reference = pd.DataFrame(columns=scdata.var_names, dtype='float32')
    label = pd.DataFrame(columns=["annotation"], dtype='str')
    for type in uni_types:
        # 实例化一个模型
        model = VAE(count_size=count_size, h_dim=h_dim, z_dim=z_dim).to(device)

        # 创建优化器
        optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)
        print("Calculate examples for {0}".format(type))
        zidx = np.where(scdata.obs['annotation'] == type)[0]
        dataset = scdataset[zidx, :]
        dataset = torch.from_numpy(np.array(dataset))
        batch_size = len(dataset)
        # 数据加载器
        data_loader = Data.DataLoader(dataset=dataset,
                                      batch_size=batch_size,
                                      shuffle=True)
        for epoch in range(num_epochs):
            kl_divs = []
            for i, x in enumerate(data_loader):
                # 获取样本，并前向传播
                x = x.to(device).view(-1, count_size)
                x_reconst, mu, log_var = model(x)

                reconst_loss = F.binary_cross_entropy(x_reconst, x, reduction='mean')
                kl_div = - 0.5 * torch.sum(1 + log_var - mu.pow(2) - log_var.exp())

                # 反向传播和优化
                loss = reconst_loss + kl_div
                optimizer.zero_grad()
                loss.backward()
                optimizer.step()

                kl_divs.append(kl_div)
                if epoch % 100 == 0 and epoch > 0:
                    print("Epoch[{}/{}] Reconst Loss: {}, KL Div: {}"
                          .format(epoch, num_epochs, reconst_loss.item(), kl_div.item()))
            if sum(kl_divs) < epoch_thresh:
                break
        print("Epoch has reach an end.")
        with torch.no_grad():
            # 重构的图像
            out, _, _ = model(x)
            npout = out.numpy()
            oidx = np.random.choice(np.arange(npout.shape[0]), size=ref_size, replace=False)
            ref = pd.DataFrame(npout[oidx, :], columns=scdata.var_names)
            # 保存参考数据集和annotation
            reference = pd.concat([reference, ref], ignore_index=True)
            types = pd.DataFrame(np.repeat(type, ref_size), columns=['annotation'])
            label = pd.concat([label, types], ignore_index=True)
    reference = scaler.inverse_transform(reference)
    reference = np.round(np.exp(reference)-1)
    reference = pd.DataFrame(reference, columns=scdata.var_names, dtype='int32')
    return reference, label

def Save_tsv_from_ref(reference: pd.DataFrame, label: pd.DataFrame, tempdir, name):
    tsvPath = tempdir + "/references/" + name + '/'
    if not os.path.exists(tsvPath):
        os.mkdir(tsvPath)
    reference.to_csv(tsvPath + '/cnt_data.tsv', sep='\t')
    label.columns = pd.Index(["bio_celltype"])
    label.to_csv(tsvPath + '/mta_data.tsv', sep='\t')
