#!/bin/bash
#
#SBATCH --job-name de
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=6GB

source /local10G/rfhorns/resources/anaconda2/bin/activate /local10G/rfhorns/resources/anaconda2/
python /local10G/rfhorns/FlyBrain/rnaseq/analysis/v7/GH146_UniqueMarkers/calc_DE_mannwhitneyu.py $1 $2 $3 $4
