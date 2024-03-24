import numpy as np
from sklearn.decomposition import PCA

from glob import glob
import sys, os

# Usage: python get_PCA_embedding.py input.a3m
# Output: df_with_PCs.json.zip

def read_msa(fil,numeric=True):
    names, seqs=[],[]
    with open(fil,'r') as f:
        for lin in f.readlines():
            if not lin.startswith('>'):
                seq = ''.join([x for x in lin if not x.islower()]).replace('\n','')
                seqs.append(seq)
            elif lin.startswith('>'):
                names.append(lin.replace('>','').replace('\n',''))

    if not numeric:
        return names, seqs
    else:
        chars='-ACDEFGHIKLMNPQRSTVWYX'
        lst=[]
        for seq in seqs:
            lst.append([chars.index(x) for x in seq])
        return names, np.array(lst)

def one_hot(x,cat=None):
    if cat is None: cat = np.max(x)+1
    oh = np.concatenate((np.eye(cat),np.zeros([1,cat])))
    return oh[x]

names, seqs = read_msa(sys.argv[1])

msa = one_hot(seqs)
msa = msa.reshape(msa.shape[0],-1)
mdl = PCA()
pca_c = mdl.fit_transform(msa)

df = pd.DataFrame({'ID': names, 'PC 1': pca_c[:,0],'PC 2': pca_c[:,1],
    'PC 3': pca_c[:,2], 'PC 4': pca_c[:,3]})

df.to_json('df_with_PCs.json.zip')

