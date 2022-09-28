import mdtraj as md
from sklearn.decomposition import PCA
import numpy as np
import os
import pandas as pd

def get_contacts(pdb_file, i=0, f=-1):
    obj = md.load_pdb(pdb_file)
    distances, pairs = md.compute_contacts(obj, scheme='ca')
    arr = md.geometry.squareform(distances, pairs)[0]
    arr = arr[i:f, i:f]
    arr = arr.flatten()
    return arr

df = pd.read_json('2cvba_feats_with_tails.json.zip')


df['contacts'] = df.apply(lambda row: get_contacts('preds_3r/'+os.path.basename(row['pdb']),i=16,f=161), axis=1)
mdl = PCA(n_components=5)
embedding = mdl.fit_transform(np.array(df.contacts.tolist()))
for i in range(5):
    df['PC %d' % (i+1)] = embedding[:,i]

df = df.drop(columns=['contacts'])
df.to_json('2cvba_feats.json.zip')

df = pd.read_json('1r26A_feats_with_tails.json.zip')

df['contacts'] = df.apply(lambda row: get_contacts('preds_3r/'+os.path.basename(row['pdb']),i=16), axis=1)
mdl = PCA(n_components=5)
embedding = mdl.fit_transform(np.array(df.contacts.tolist()))
for i in range(5):
    df['PC %d' % (i+1)] = embedding[:,i]

df = df.drop(columns=['contacts'])
df.to_json('1r26a_feats.json.zip')