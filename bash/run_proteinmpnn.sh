#!/bin/bash
#SBATCH --job-name=5F91_seq
#SBATCH --cpus-per-task=8
#SBATCH --time=100:00:00
#SBATCH --mem=100g
#SBATCH --gres=gpu:p100:1,lscratch:100
#SBATCH -p gpu

# If you are running on HPC and the HPC has already compiled ProteinMPNN
# module purge
# module load ProteinMPNN
conda activate proteinMPNN

folder_with_pdbs="./results/rfdiffusion_pdb"
output_dir="./results/proteinmpnn_sequences"

# the number of sequences to design for each model
NUMBER_SEQ_PER_TARGET=10

mkdir -p "$output_dir"

path_for_parsed_chains="${output_dir}/parsed_pdbs.jsonl"
path_for_assigned_chains="${output_dir}/assigned_pdbs.jsonl"

CHAIN_TO_DESIGN="B"   # peptide chain only

python ./ProteinMPNN/helper_scripts/parse_multiple_chains.py \
  --input_path="$folder_with_pdbs" \
  --output_path="$path_for_parsed_chains"

python ./ProteinMPNN/helper_scripts/assign_fixed_chains.py \
  --input_path="${path_for_parsed_chains}" \
  --output_path="${path_for_assigned_chains}" \
  --chain_list "${CHAIN_TO_DESIGN}"

python ./ProteinMPNN/protein_mpnn_run.py \
  --jsonl_path "$path_for_parsed_chains" \
  --chain_id_jsonl "$path_for_assigned_chains" \
  --out_folder "$output_dir" \
  --num_seq_per_target ${NUMBER_SEQ_PER_TARGET} \
  --sampling_temp "0.1" \
  --seed 37 \
  --batch_size 1
