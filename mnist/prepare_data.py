import numpy as np
import torch as t
import torchvision as tv
from openTSNE import TSNE
from sklearn.decomposition import PCA

from config import raw_path, labels_path, pca_path, tsne_path, N


def load_tsne():
    if tsne_path.exists():
        print("==== Loading tSNE Embedding ====\n")
        proj_tsne = np.load(tsne_path)
    else:
        proj_pca = load_pca()
        print("==== Fitting tSNE Embedding ====")
        tsne = TSNE(
            n_jobs=4,
            random_state=20260417,
            verbose=True,
        )
    
        proj_tsne = tsne.fit(proj_pca)
        np.save(tsne_path, proj_tsne)

    return proj_tsne

def load_pca():
    if pca_path.exists():
        print("==== Loading PCA Embedding ====")
        proj_pca = np.load(pca_path)
    else:
        sampleX_t = load_raw()
        print("==== Fitting PCA Embedding ====")
        pca = PCA(n_components=100)
        proj_pca = pca.fit_transform(sampleX_t.reshape(-1,784))
        np.save(pca_path, proj_pca)

    return proj_pca

def load_raw():
    if raw_path.exists():
        print("==== Loading Raw Data ====")
        sampleX_t = np.load(raw_path)
        sampleY_t = np.load(labels_path)
    else:
        print("==== Preparing Raw Data ====")
        ## Set up Data tensor
        transform = tv.transforms.Compose([tv.transforms.ToTensor()])
        
        train_data = tv.datasets.MNIST("data/train",train=True, download=True, transform = transform)
        test_data = tv.datasets.MNIST("data/test",train=False, download=True, transform=transform)
        
        kwargs = {'num_workers': 1, 'pin_memory': False}
        
        train_loader = t.utils.data.DataLoader(train_data,batch_size=100,shuffle=False,**kwargs)
        test_loader = t.utils.data.DataLoader(test_data,batch_size=100,shuffle=False,**kwargs)
        
        sampleX = []
        sampleY = []
        for loader in [train_loader, test_loader]:
            for b in loader:
                bx, by = b
                sampleX.append(bx)
                sampleY.append(by)
    
        sampleX_t = np.array(sampleX).reshape([N,1,28,28])
        sampleY_t = np.array(sampleY).reshape([N])
        np.save(labels_path, sampleY_t)    
        np.save(raw_path, sampleX_t)

    return sampleX_t, sampleY_t


def load_labels():
    if raw_path.exists():
        sampleY_t = np.load(labels_path)
    else:
        _, sampleY_t = load_raw()
    
    return sampleY_t
