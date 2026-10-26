import numpy as np
import torch as t
import torchvision as tv
from openTSNE import TSNE
from sklearn.decomposition import PCA

from config import raw_path, labels_path, pca_path, tsne_path, N

## Pull out raw data and prepare as Numpy array
if raw_path.exists():
    print("Loading Raw")
    sampleX_t = np.load(raw_path)
    sampleY_t = np.load(labels_path)
else:
    print("Preparing raw")
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

## Perform PCA
if pca_path.exists():
    print("Loading PCA")
    proj_pca = np.load(pca_path)
else:
    print("Fitting PCA")
    pca = PCA(n_components=100)
    proj_pca = pca.fit_transform(sampleX_t.reshape(-1,784))
    np.save(pca_path, proj_pca)


## Perform tSNE
if tsne_path.exists():
    print("Loading tsne")
    proj_tsne = np.load(tsne_path)
else:
    print("Fitting tsne")
    tsne = TSNE(
        n_jobs=4,
        random_state=20260417,
        verbose=True,
    )

    proj_tsne = tsne.fit(proj_pca)
    np.save(tsne_path, proj_tsne)
