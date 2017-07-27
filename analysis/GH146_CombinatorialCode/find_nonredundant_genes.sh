#!/bin/bash
#
#SBATCH --job-name find_nonredundant_genes
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=6GB

source /local10G/rfhorns/resources/anaconda2/bin/activate /local10G/rfhorns/resources/anaconda2/
python /local10G/rfhorns/FlyBrain/rnaseq/analysis/v7/GH146_CombinatorialCode_RevisedGeneLists/find_nonredundant_genes.py
