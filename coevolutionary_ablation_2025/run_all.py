import numpy as np
#import string
import pandas as pd
import sys, os, re

from colabdesign import mk_af_model, clear_mem
from colabdesign.af.contrib import predict
from colabdesign.shared.protein import _np_rmsd

from utils import run

autoinb=['049','056','021','014','045']
active=['001','002','005','024','025']
rest=['000', '007', '011', '012', '015', '016', '017', '019', '022', '028', '032', '035', '038', '039', '041', '050', '054', '057', '059', '060', '066', '083', '085', '091', '093', '096', '099', '108', '120', '124', '127', '138', '143', '147', '151', '157', '173', '194', '238']
# more than 10 seqs not including query and not in autoinb or active

sequence = "MQSWYLLYCKRGQLQRAQEHLERQAVNCLAPMITLEKIVRGKRTAVSEPLFPNYLFVEFDPEVIHTTTINATRGVSHFVRFGASPAIVPSAVIHQLSVYKPKDIVDPATPYPGDKVIITEGAFEGFQAIFTEPDGEARSMLLLNLINKEIKHSVKNTEFRKA" 

jobname = "RfaH_autoinb"
run(sequence,jobname,'RfaH_msas/RFAH',autoinb,scramble_cols=True,verbose=True)
run(sequence,jobname,'RfaH_msas/RFAH',autoinb,verbose=True)

jobname = "RfaH_active"
run(sequence,jobname,'RfaH_msas/RFAH',active,scramble_cols=True,verbose=True)
run(sequence,jobname,'RfaH_msas/RFAH',active,verbose=True)

jobname = "RfaH_rest"
run(sequence,jobname,'RfaH_msas/RFAH',rest,scramble_cols=True,verbose=True,num_seeds=4)
run(sequence,jobname,'RfaH_msas/RFAH',rest,verbose=True,num_seeds=4)

jobname = "RfaH_rest_nodropout"
run(sequence,jobname,'RfaH_msas/RFAH',rest,scramble_cols=True,verbose=True,num_seeds=4,use_dropout=False)
run(sequence,jobname,'RfaH_msas/RFAH',rest,verbose=True,num_seeds=4,use_dropout=False)

s2h=['000','046','008']
duj=['020','055','083']
rest=['001', '002', '004', '007', '009', '013', '017', '021', '024', '027', '028', '030', '032', '038', '039', '040', '041', '042', '044', '052', '060', '064', '069', '073', '074', '078']

msa_path='Mad2_msas/1S2H'

sequence = "GMALQLSREQGITLRGSAEIVAEFFSFGINSILYQRGIYPSETFTRVQKYGLTLLVTTDLELIKYLNNVVEQLKDWLYKCSVQKLVVVISNIESGEVLERWQFDIECDKTAKDDSAPREKSQKAIQDEIRSVIAQITATVTFLPLLEVSCSFDLLIYTDKDLVVPEKWEESGPQFITNSEEVRLRSFTTTIHKVNSMVAYKIPVND" 

jobname = "Mad2_1S2H"
run(sequence,jobname,msa_path,s2h,scramble_cols=True,verbose=True)
run(sequence,jobname,msa_path,s2h,verbose=True)

jobname = "Mad2_1DUJ"
run(sequence,jobname,msa_path,duj,scramble_cols=True,verbose=True)
run(sequence,jobname,msa_path,duj,verbose=True)

jobname = "Mad2_rest"
run(sequence,jobname,msa_path,rest,scramble_cols=True,verbose=True,num_seeds=4)
run(sequence,jobname,msa_path,rest,verbose=True,num_seeds=4)
