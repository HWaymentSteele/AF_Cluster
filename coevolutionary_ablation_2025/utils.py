
import numpy as np
import pandas as pd
import sys, os, re

from colabdesign import mk_af_model, clear_mem
from colabdesign.af.contrib import predict
from colabdesign.shared.protein import _np_rmsd

def scramble_msa_cols(msa, dtx, seed=0):
    N, L = msa.shape

    msa_to_scramble = msa[1:]
    dtx_to_scramble = dtx[1:]
    new_cols_msa, new_cols_dtx=[],[]
    
    for pos in range(L):
        np.random.seed(seed=pos+seed)
        order = np.random.choice(range(N-1),N-1)
        new_cols_msa.append(msa_to_scramble[order,pos])
        new_cols_dtx.append(dtx_to_scramble[order,pos])

    new_msa = np.hstack([msa[0].reshape(-1,1), np.vstack(new_cols_msa)]).T
    new_dtx = np.hstack([dtx[0].reshape(-1,1), np.vstack(new_cols_dtx)]).T

    return new_msa, new_dtx

def run(sequence, jobname, msa_path, msa_cluster_list, scramble_cols=False, verbose=False, num_seeds=4,use_dropout=False):
    copies = 1 
    msa_method = "mmseqs2"
    pair_mode = "unpaired_paired" 
    
    qid = 0 
    do_not_filter = False 
    min_samples_in_cluster = 10 
    
    template_mode = "none" 
    use_templates = template_mode in ["mmseqs2","custom"]
    pdb = "" 
    chain = "" 
    flexible = False 
    propagate_to_copies = True
    rm_interchain = False
    rm_sidechain = rm_sequence = flexible
    
    # process sequence
    sequences = sequence.split(":")
    u_sequences = predict.get_unique_sequences(sequences)
    if len(sequences) > len(u_sequences):
        print("WARNING: use copies to define homooligomers")
    u_lengths = [len(s) for s in u_sequences]
    
    sub_seq = "".join(u_sequences)
    seq = sub_seq * copies
    
    jobname = f"{jobname}_{predict.get_hash(seq)[:5]}"
    input_opts = {"sequence":u_sequences,
                  "copies":copies,
                  "msa_method":msa_method,
                  "pair_mode":pair_mode,
                  "do_not_filter":do_not_filter,
                  "template_mode":template_mode,
                  "propagate_to_copies":propagate_to_copies}
    
    os.makedirs(jobname, exist_ok=True)
    clustered_msas, clustered_dtxs = [],[]

    for s in msa_cluster_list:
        x=f'{msa_path}_{s}.a3m'
        m, d = predict.parse_a3m(x)

        clustered_msas.append(m)
        clustered_dtxs.append(d)
        
    batches = [None]
    model_type = "auto" 
    model = "all"
    num_recycles = 3
    recycle_early_stop_tolerance = 0.0
    use_initial_guess = False
    num_msa = 512 
    num_extra_msa = 1024 
    use_cluster_profile = False 
    use_mlm = False 
    use_dropout = use_dropout 
    
    if scramble_cols:
        scramble_print='_SCRAMBLE'
    else:
        scramble_print=''
    seed = 0 
    
    if model_type == "monomer (ptm)":
        use_multimer = False
    elif model_type == "multimer (v3)":
        use_multimer = True
    elif len(u_lengths) > 1 or copies > 1:
        use_multimer = True
    else:
        use_multimer = False
    
    model_opts = {"num_msa":num_msa, # number of sequences to use
                  "num_extra_msa":num_extra_msa,
                  "num_templates":len(batches),
                  "use_mlm":True,
                  "use_cluster_profile":use_cluster_profile,
                  "use_multimer":use_multimer,
                  "use_templates":use_templates,
                  "use_batch_as_template":False,
                  "use_dgram":True,
                  "protocol":"hallucination",
                  "best_metric":"pae",
                  "optimize_seq":False,
                  "debug":False,
                  "clear_prev":False,
                 "data_dir":'/ssd/NMR_project_data/'}
    
    if "af" in dir():
      if model_opts != model_opts_:
        if model_opts["use_multimer"] == af._args["use_multimer"] \
        or model_opts["use_templates"] == af._args["use_templates"]:
          old_params = dict(zip(af._model_names,
                                af._model_params))
        else:
          print("loading alphafold params")
          old_params = {}
        af = mk_af_model(old_params=old_params,
                         **model_opts)
        model_opts_ = predict.copy_dict(model_opts)
    else:
      print("loading alphafold params")
      af = mk_af_model(**model_opts)
      model_opts_ = predict.copy_dict(model_opts)
    
    run_opts = {"seed":seed,
                "use_mlm":use_mlm,
                "use_dropout":use_dropout,
                "num_recycles":num_recycles,
                "model":model,
                "use_initial_guess":use_initial_guess}
    
    af.prep_inputs(u_lengths, copies=copies, seed=seed)
    if use_templates:
      af.set_opt(use_initial_guess=use_initial_guess)
      for n,batch in enumerate(batches):
        af.set_template(batch=batch, n=n)
      af.set_opt("template",
                 rm_sc=rm_sidechain,
                 rm_seq=rm_sequence,
                 rm_ic=rm_interchain)
    af.set_opt("mlm",
               replace_fraction=0.15 if use_mlm else 0.0)
    
    if model == "all":
      models = af._model_names
    else:
      models = [af._model_names[int(model) - 1]]
    
    pdb_path = f"{jobname}/pdb"
    os.makedirs(pdb_path, exist_ok=True)
    
    seeds = list(range(seed,seed+num_seeds))
    print("running prediction")
    
    outputs=[]
    
    for c_ind in range(len(clustered_msas)):
    
      for seed in seeds:
        af.set_seed(seed)
        if scramble_cols:
            local_msa_ , local_dtx_ = scramble_msa_cols(clustered_msas[c_ind], clustered_dtxs[c_ind],seed=seed)
            af.set_msa(local_msa_ , local_dtx_)
    
        else:
          af.set_msa(clustered_msas[c_ind], clustered_dtxs[c_ind])
            
        for model in models:
          o={}
          recycle = 0
          af._inputs.pop("prev",None)
          stop_recycle = False
          while recycle < num_recycles + 1:
            af.predict(dropout=use_dropout, models=[model], verbose=False)
    
            print_str = f"cluster={c_ind} seed={seed} model={model} recycle={recycle}"
            print_key = ["plddt","ptm"]
            if len(af._lengths) > 1: print_key.append("i_ptm")
            for k in print_key:
              print_str += f" {k}={af.aux['log'][k]:.3f}"
    
            af._inputs["prev"] = af.aux["prev"]
            af._save_results(save_best=True, verbose=False)
            af._k += 1
            output_pdb = f"{pdb_path}/{model}_cluster{c_ind}_r{recycle}_seed{seed}{scramble_print}.pdb"
            af.save_current_pdb(output_pdb)
    
            recycle += 1
            if recycle > 1:
              rmsd_tol = _np_rmsd(af._tmp["traj"]["xyz"][-2],
                                  af._tmp["traj"]["xyz"][-1],
                                  use_jax=False)
              if rmsd_tol < recycle_early_stop_tolerance:
                stop_recycle = True
              print_str += f" rmsd_tol={rmsd_tol:.3f}"
            if verbose: print(print_str)
            if stop_recycle: break
            for k in print_key:
              o.update({k: af.aux['log'][k]})
    
            o.update({'pdb_path': output_pdb})
            o.update({'model': model})
            o.update({'cluster': c_ind})
            o.update({'cluster_size': clustered_msas[c_ind].shape[0]-1})
            outputs.append(o)
    
    outputs = pd.DataFrame.from_records(outputs)
    outputs.to_csv(f'{jobname}{scramble_print}.csv')
