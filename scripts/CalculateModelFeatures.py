from glob import glob
import sys, os, argparse
import pandas as pd
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
from polyleven import levenshtein

from pandarallel import pandarallel

pandarallel.initialize(use_memory_fs=False)
lst=[]

def read_b_factor(pdb_file):
    vals=[]

    with open(pdb_file,'r') as f:
        for lin in f.readlines()[1:-3]:
            fields = lin.split()
            if fields[2] =='CA':
                vals.append(float(fields[10]))
    return vals

def calc_rmsd(decoy, ref_obj, debug=False):
    if isinstance(decoy,str):
        decoy = md.load_pdb(decoy)

    ref_atoms = [x.index for x in ref_obj.top.residues]
    ref_seq = ref_obj.top.to_fasta()[0]
    decoy_seq = decoy.top.to_fasta()[0]

    string='name CA and ('+ ' or '.join("resid %d" % x for x in ref_atoms)+')'

    if debug:
        
        print(ref_seq)
        print(decoy_seq)
        if len(decoy.top.select(string)) != len(ref_obj.top.select(string)):
            print(len(decoy.top.select(string)), len(ref_obj.top.select(string)), string)

        print([x.index for x in ref_obj.top.residues])
        print([x.index for x in decoy.top.residues])

        print(ref_obj.top.select(string))
        print(decoy.top.select(string))
    return md.rmsd(decoy, ref_obj,  atom_indices=decoy.top.select(string), ref_atom_indices=ref_obj.top.select(string))[0]

def calc_dssp(decoy):
    if isinstance(decoy,str):
        decoy = md.load_pdb(decoy)

    return ''.join([x for x in md.compute_dssp(decoy).tolist()[0]])


if __name__=='__main__':
    p = argparse.ArgumentParser(description=
    """
    Calculate structure features of decoys.

    H Wayment-Steele, 2022
    """)

    p.add_argument("input_pdbs", nargs='*', action='store',help='Path to pdbs to featurize.')
    p.add_argument("-o", action="store", default='pdb_features.json.zip', help='name of output file.')
    p.add_argument("--test", action='store_true', help='Tests first 10 decoys.')
    p.add_argument('--ref_struct', action='store', help='PDB reference structure.')
    p.add_argument('--plot_results', action='store', help='Plot results.')

    args = p.parse_args()

    lst=[]
    print("Found %d pdbs" % len(args.input_pdbs))

    df = pd.DataFrame({'pdb': args.input_pdbs})

    if args.test:
        df = df.iloc[:10]


    df['pLDDT_vector'] = df.parallel_apply(lambda row: read_b_factor(row['pdb']), axis=1)
    df['mean_pLDDT'] = df.parallel_apply(lambda row: np.mean(row['pLDDT_vector']), axis=1)
    df['dssp_string'] = df.parallel_apply(lambda row: calc_dssp(row['pdb']), axis=1)

    if args.ref_struct is not None:
        pdb_name = os.path.basename(args.ref_struct).replace('.pdb','')
        df['rmsd_ref'] = df.parallel_apply(lambda row: calc_rmsd(row['pdb'], ref_obj), axis=1)

        ref_obj = md.load_pdb(args.ref_struct)
        ref_dssp = calc_dssp(ref_obj)
        df['dist_dssp'] = df.parallel_apply(lambda row: levenshtein(row['dssp_string'], ref_dssp), axis=1)
        
    df.head()

    if args.plot_results and args.ref_struct is not None:
        plt.figure(figsize=(8,4))
        plt.subplot(1,2,1)
        plt.scatter(df['rmsd_ref'].to_numpy(), df['mean_pLDDT'].to_numpy())
        for _, row in df.iterrows():
            plt.text(row['rmsd_ref'], row['mean_pLDDT'], os.path.basename(row['pdb'].replace('.pdb','').split('_')[-1]))
        plt.xlabel('RMSD to ref')
        plt.ylim([0,100])
        plt.title(pdb_name)
        plt.ylabel('pLDDT')

        plt.subplot(1,2,2)
        plt.scatter(df['dist_dssp'].to_numpy(), df['mean_pLDDT'].to_numpy())
        for _, row in df.iterrows():
            plt.text(row['dist_dssp'], row['mean_pLDDT'], os.path.basename(row['pdb'].replace('.pdb','').split('_')[-1]))
        plt.xlabel('dist plddt to ref')
        plt.ylim([0,100])
        plt.xlim([0, len(ref_dssp)])

        plt.ylabel('pLDDT')
        plt.savefig('%s_plddt_vs_rmsd.pdf' % pdb_name,bbox_inches='tight')

    df.to_json(args.o)
    print('Wrote to', args.o)

