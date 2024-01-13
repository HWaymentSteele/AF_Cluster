# AF2 commands to replicate structure models in AF-Cluster paper

Hannah Wayment-Steele, last updated Jan 2, 2024

For experiments performed in ColabFold, path to zipped files of ColabFold output are provided.
The `run_af2.py` script, which all the non-ColabFold predictions used, is contained in this repo at `/scripts/run_AF2.py` .

## Fig 1c

Left: default KaiBTE MSA in ColabFold, 3 recycles, depicting highest pLDDT of 5 models. 
```data_sep2022/colabfold_outputs/2qke_orig_default_colabfold_jul2022.result.zip```

Right: Closest 50 sequences of KaiBTE MSA in ColabFold, 3 recycles, depicting highest pLDDT of 5 models.  
```data_sep2022/colabfold_outputs/2qke_top50_seqs_jul2022.result.zip```

## Fig 1d,e,f

AF-Cluster models of KaiBTE. Generated using run_af2.py, model 1, 3 recycles, MSAs from AF-Cluster routine (`ClusterMSA.py`)

```python run_af2.py AF_Cluster/data_sep2022/AF_Cluster/data_sep2022/00_KaiB/msas/2QKEE_00.a3m --model_num 1 --recycles 3 --output_dir .```

## Fig 2b

Models of 4 KaiB variants. Generated using run_af2.py, model 1, 3 recycles, “closest-10” MSAs (closest 10 sequences in phylogenetic tree by edit distance).

Note that this MSA was not generated in the ColabFold/MMseqs routine, but rather in the process of constructing the KaiB phylogenetic tree, detailed in the paper Methods under "Phylogenetic tree construction".

```python run_af2.py AF_Cluster/data_sep2022/02_KaiB_fromPhylTree/msas/WP_069333443.1.a3m --model_num 1 --recycles 3 --output_dir .```

R. sphaeroides:	WP_069333443

T. elongatus:	WP_011056333

S. elongatus:	WP_011242647

L. pneumophilia:	WP_027265801

## Fig 2c,e

KaiBTV-4 model, generated using ColabFold and MSA of 10 closest sequences in phylogenetic tree. Model 1, 12 recycles.
```data_sep2022/colabfold_outputs/kaibtv4_local_MSA.result.zip```

## Fig 3b

Models of point mutants: generated using single sequences, model=1, 12 recycles

Fastas of single sequences are in pt_mutant_fastas_apr7.zip. First unzip:
```unzip pt_mutant_fastas_apr7.zip```

```python run_af2.py AF_Cluster/data_sep2022/03_KaiBRS_and_pt_muts/indiv_fastas_apr7/V83D_Y35C.fasta --model_num 1 --recycles 12 --output_dir .```

## Fig 4b

RfaH models. Generated using `run_af2.py`, model 1, 3 recycles, MSAs from AF-Cluster routine (`/scripts/ClusterMSA.py`)

```python run_af2.py AF_Cluster/data_sep2022/04_OtherFoldswitchers/00_RfaH/msas/RFAH_000.a3m --model_num 1 --recycles 3 --output_dir .```

## Fig 4d

Mad2 models. Generated using `run_af2.py`, model 1, 3 recycles, MSAs from AF-Cluster routine (`/scripts/ClusterMSA.py`)

```python run_af2.py AF_Cluster/data_sep2022/01_Mad2/1S2H_msas/1S2H_000.a3m --model_num 1 --recycles 3 --output_dir .```

## Fig 5b

AF-Cluster initial screen

Generated using `run_af2.py`, model 1, 1 recycle, MSAs from AF-Cluster routine (ClusterMSA.py)

```python run_af2.py AF_Cluster/data_sep2022/05_Thioredoxins/00_Mpt53/msas/1LU4A_0.a3m --model_num 1 --recycles 1 --output_dir .```

## Fig 5c

More sampling on Mpt53

Generated using `run_af2.py`, model 1, 3 recycles, MSAs from AF-Cluster routine (ClusterMSA.py)

```python run_af2.py AF_Cluster/data_sep2022/05_Thioredoxins/00_Mpt53/msas/1LU4A_0.a3m --model_num 1 --recycles 3 --output_dir .```

## Extended Data Fig 1a 

KaiB-TE (PDB: 2QKE), closest 50, closest 50-100

Generated using ColabFold, 3 recycles. 5 models depicted.

Closest 50: ```data_sep2022/colabfold_outputs/ColabFold output: data_sep2022/colabfold_outputs/2qke_top50_seqs_jul2022.result.zip```

Closest 50-100: ```data_sep2022/colabfold_outputs/2qke_next50_aug2023.result.zip```

## Extended data Fig 2d, 2e

Sweep of eps values, same method as Fig 1d,e,f

## Extended Data Fig 3d

KaiB AF-Cluster models: same method/data as Fig 1d,e,f

## Extended Data Fig 5a,b

Mutant models: same method/data as Fig 3b

## Extended data Fig 6a

RfaH in AF2, complete MSA

Generated using ColabFold, 3 recycles, depicting 5 models.

```data_sep2022/colabfold_outputs/Rfah_3_recycles_default_colabfold.result.zip```

## Extended data fig 6c

Other fold-switchers

Selecase:

```python run_af2.py AF_Cluster/data_sep2022/02_selecase/msas/4QHF_000.a3m --model_num 1 --recycles 3 --output_dir .```

Lymphotactin:

```python run_af2.py AF_Cluster/data_sep2022/03_lymphotactin/2JP1_msas/2JP1_000.a3m --model_num 1 --recycles 3 --output_dir .```

CLIC1:

```python run_af2.py AF_Cluster/data_sep2022/04_OtherFoldswitchers/04_CLIC1/1K0N_000.a3m --model_num 1 --recycles 3 --output_dir .```

## Extended data Fig 7a,b

Enhanced sampling Mpt53: Same method as Fig. 5c

## Extended data Fig 8

Homologues to Mpt53: Same method as Fig. 5c

## Extended data Fig 10c

G_A/G_B homologues: Performed using ColabDesign code, 0 recycles, 8 seeds, 5 models. 

Colab notebook: 
```data_sep2022/colabfold_outputs/GAGB_mutations_in_colabdesign.ipynb```

## Supplemental discussion Fig 1a,b,d

Comparing models from single sequence and “closest-10” MSAs.

“Closest-10” models are from the methods described in Fig. 2b (3 recycles, model 1).
Single-sequence models were generated analogously (3 recycles, model 1).

First unzip single seq fastas:
```unzip AF_Cluster/data_sep2022/00_KaiB/fastas_single_seq.zip```

```python run_af2.py  AF_Cluster/data_sep2022/00_KaiB/fastas_single_seq/WP_069333443.1.fasta --model_num 1 --recycles 3 --output_dir .```

## Supplemental discussion Fig 1c

10 tests of shuffling residues within columns.
These were performed in ColabFold with 3 recycles and model 1. 
```data_sep2022/colabfold_outputs/Supp_Discussion_fig1c_scrambled_msas/*.result.zip```

Supplemental figure 2: Testing increased sampling on KaiBRS.
Performed in colabfold, using 16 seeds, 3 recycles, 5 models, use_dropout=True.
```data_sep2022/colabfold_outputs/kaibrs_lots_of_sampling_supp_discussion_fig2.result.zip```

Note: in this .zip, .pickle files were not uploaded due to file size.

## Followup on KaiBTV-4

Testking KaiBTV-4 predictions across all models and recycles up to 12. Below, \<NUM\> = 1-5, \<N_RECYCLES\>=0-12.

"local 10" MSA:
```python run_af2.py data_sep2022/06_kaibtv4_followup/msas/kaibtv4_local10.a3m --model_num <NUM> --recycles <N_RECYCLES> --output_dir .```

Single sequence mode:
```python run_af2.py data_sep2022/06_kaibtv4_followup/msas/kaibtv4_single_seq.a3m --model_num <NUM> --recycles <N_RECYCLES> --output_dir .```

Column-scrambled replicate, \<REP\>=0-9:
```python run_af2.py data_sep2022/06_kaibtv4_followup/msas/scrambled/scramble_<REP>.a3m --model_num <NUM> --recycles <N_RECYCLES> --output_dir .```
