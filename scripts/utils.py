import numpy as np
from Bio import SeqIO

def lprint(string,f):
    print(string)
    f.write(string+'\n')

def load_fasta(fil):
    seqs, IDs =[], []
    with open(fil) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                seq = ''.join([x for x in record.seq])
                IDs.append(record.id)
                seqs.append(seq)
    return IDs, seqs

def write_fasta(names, seqs, outfile='tmp.fasta'):
        with open(outfile,'w') as f:
                for nm, seq in list(zip(names, seqs)):
                        f.write(">%s\n%s\n" % (nm, seq))

def dihedral_wrapper(traj):
    """Featurize an MD trajectory into a vector space via calculation
    of dihedral (torsion) angles of alpha carbon backbone
    
    Lifted from MSMBuilder, RT McGibbon

    Parameters
    ----------
    traj : mdtraj.Trajectory
        A molecular dynamics trajectory to featurize.

    Returns
    -------
    features : np.ndarray, dtype=float, shape=(n_samples, n_features)
        A featurized trajectory is a 2D array of shape
        `(length_of_trajectory x n_features)` where each `features[i]`
        vector is computed by applying the featurization function
        to the `i`th snapshot of the input trajectory.

    """

    ca = [a.index for a in traj.top.atoms if a.name == 'CA']
    if len(ca) < 4:
        return np.zeros((len(traj), 0), dtype=np.float32)

    alpha_indices = np.array(
        [(ca[i - 1], ca[i], ca[i+1], ca[i + 2]) for i in range(1, len(ca) - 2)])
    result = md.compute_dihedrals(traj, alpha_indices)

    return result[0]
    
def consensusVoting(seqs):
    ## Find the consensus sequence
    consensus = ""
    residues = "ACDEFGHIKLMNPQRSTVWY-"
    n_chars = len(seqs[0])
    for i in range(n_chars):
        baseArray = [x[i] for x in seqs]
        baseCount = np.array([baseArray.count(a) for a in list(residues)])
        vote = np.argmax(baseCount)
        consensus += residues[vote]

    return consensus

def encode_seqs(seqs, max_len=108, alphabet=None):
    
    if alphabet is None:
        alphabet = "ACDEFGHIKLMNPQRSTVWY-"
    
    arr = np.zeros([len(seqs), max_len, len(alphabet)])
    for j, seq in enumerate(seqs):
        for i,char in enumerate(seq):
            for k, res in enumerate(alphabet):
                if char==res:
                    arr[j,i,k]+=1
    return arr.reshape([len(seqs), max_len*len(alphabet)])
