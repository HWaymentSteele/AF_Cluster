#!/bin/bash
sbatch <<EOT
#!/bin/bash

#SBATCH --partition=gpu,gpu_requeue
#SBATCH --gres=gpu:1
#SBATCH --time=2:00:00
#SBATCH --mem=4G
#SBATCH --output=$PWD/slurm_out/slurm-%j.out
#SBATCH --job-name=${1}
###
source ~/.bash_profile
conda activate AF
module load cuda

export XLA_FLAGS=--xla_gpu_cuda_data_dir=/n/sw/helmod-rocky8/apps/Core/cuda/12.2.0-fasrc01/cuda
cd $PWD

python run_af2.py ${1} --model_num ${2} --recycles ${3} --output_dir $PWD/preds/
EOT
# first arg: local dir containing files to run
# second arg: model number to run
# third arg: n recycles
