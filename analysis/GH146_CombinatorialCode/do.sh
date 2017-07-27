NUM_REPLICATES=100

for i in $(seq 1 $NUM_REPLICATES)
do
    sbatch find_nonredundant_genes.sh
done
