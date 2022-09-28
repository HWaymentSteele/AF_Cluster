import sys
import os
import argparse
import hashlib
import jax
import jax.numpy as jnp
import numpy as np
import re
import subprocess
from glob import glob
import pandas as pd
from collections import namedtuple

parser = argparse.ArgumentParser()
parser.add_argument("input_msa", help="MSA to run prediction. First seq will be the predicted seq")
parser.add_argument("--recycles", type=int, default=1, help="Number of recycles when predicting")
parser.add_argument("--model_num", type=int, default=1, help="Which AF2 model to use")
parser.add_argument("--seed", type=int, default=0, help="RNG Seed")
parser.add_argument("--verbose", action='store_true', help="print extra")
parser.add_argument("--deterministic", action='store_true', help="make all data processing deterministic (no masking, etc.)")
parser.add_argument("--af2_dir", default="/home/haw053/software/alphafold-2.2.0/", help="AlphaFold code and weights directory")
parser.add_argument("--output_dir", default="/home/haw053/metamorph_benchmark/", help="Where to write output files")

args = parser.parse_args()

sys.path.append(args.af2_dir)

from alphafold.model import model
from alphafold.model import config
from alphafold.model import data

from alphafold.data import parsers
from alphafold.data import pipeline

from alphafold.common import protein
from alphafold.common import residue_constants

"""
Create an AlphaFold model runner
name -- The name of the model to get the parameters from. Options: model_[1-5]
"""
def make_model_runner(name, recycles):
  cfg = config.model_config(name)      

  cfg.data.common.num_recycle = recycles
  cfg.model.num_recycle = recycles
  cfg.data.eval.num_ensemble = 1
  if args.deterministic:
    cfg.data.eval.masked_msa_replace_fraction = 0.0
    cfg.model.global_config.deterministic = True

  print('to here')
  params = data.get_model_haiku_params(name, args.af2_dir + 'data/')
  print('read params')

  return model.RunModel(cfg, params)

"""
Create a feature dictionary for input to AlphaFold
runner - The model runner being invoked. Returned from `make_model_runner`
sequence - The target sequence being predicted
templates - The template features being added to the inputs
seed - The random seed being used for data processing
"""
def make_processed_feature_dict(runner, a3m_file, name="test", templates=None, seed=0):
  feature_dict = {}

  # assuming sequence is first entry in msa

  with open(a3m_file,'r') as msa_fil:
    sequence = msa_fil.read().splitlines()[1].strip()

  feature_dict.update(pipeline.make_sequence_features(sequence, name, len(sequence)))

  with open(a3m_file,'r') as msa_fil:
    msa = pipeline.parsers.parse_a3m(msa_fil.read())

  feature_dict.update(pipeline.make_msa_features([msa]))

  if templates is not None:
    feature_dict.update(templates)
  else:
    feature_dict.update(empty_placeholder_template_features(num_templates=0, num_res=len(sequence)))


  processed_feature_dict = runner.process_features(feature_dict, random_seed=seed)

  return processed_feature_dict

"""
Make a set of empty features for no-template evalurations
"""
def empty_placeholder_template_features(num_templates, num_res):
  return {
      'template_aatype': np.zeros(
          (num_templates, num_res,
           len(residue_constants.restypes_with_x_and_gap)), dtype=np.float32),
      'template_all_atom_masks': np.zeros(
          (num_templates, num_res, residue_constants.atom_type_num),
          dtype=np.float32),
      'template_all_atom_positions': np.zeros(
          (num_templates, num_res, residue_constants.atom_type_num, 3),
          dtype=np.float32),
      'template_domain_names': np.zeros([num_templates], dtype=object),
      'template_sequence': np.zeros([num_templates], dtype=object),
      'template_sum_probs': np.zeros([num_templates], dtype=np.float32),
  }

"""
Package AlphaFold's output into an easy-to-use dictionary
prediction_result - output from running AlphaFold on an input dictionary
processed_feature_dict -- The dictionary passed to AlphaFold as input. Returned by `make_processed_feature_dict`.
"""
def parse_results(prediction_result, processed_feature_dict):
  b_factors = prediction_result['plddt'][:,None] * prediction_result['structure_module']['final_atom_mask']  
  dist_bins = jax.numpy.append(0,prediction_result["distogram"]["bin_edges"])
  dist_mtx = dist_bins[prediction_result["distogram"]["logits"].argmax(-1)]
  contact_mtx = jax.nn.softmax(prediction_result["distogram"]["logits"])[:,:,dist_bins < 8].sum(-1)

  out = {"unrelaxed_protein": protein.from_prediction(processed_feature_dict, prediction_result, b_factors=b_factors),
        "plddt": prediction_result['plddt'],
        "pLDDT": prediction_result['plddt'].mean(),
        "dists": dist_mtx,
        "adj": contact_mtx}

  out.update({"pae": prediction_result['predicted_aligned_error'],
              "pTMscore": prediction_result['ptm']})
  return out

def write_results(result, pdb_out_path):
  plddt = float(result['pLDDT'])
  ptm = float(result["pTMscore"])
  print('plddt: %.3f' % plddt)
  print('ptm  : %.3f' % ptm)

  pdb_lines = protein.to_pdb(result["unrelaxed_protein"])
  with open(pdb_out_path, 'w') as f:
    f.write(pdb_lines)

model_name = "model_{}_ptm".format(args.model_num)
results_key = model_name + "_seed_{}".format(args.seed)

runner = make_model_runner(model_name, args.recycles)

print('made model runner')

name=os.path.basename(args.input_msa).replace('.a3m','')

features = make_processed_feature_dict(runner, args.input_msa, name=name, seed=args.seed)
result = parse_results(runner.predict(features, random_seed=args.seed), features)
write_results(result, args.output_dir+name+'.pdb')



