DF_EXPR=/local10G/rfhorns/FlyBrain/rnaseq/analysis/v7/df_Mz19DevelopmentWithDL3.csv
NUM_TO_SAMPLE=12
NUM_REPLICATES=100

# DL3 vs DA1

NAMES1=/local10G/rfhorns/FlyBrain/rnaseq/analysis/v7/names_DL3.txt
NAMES2=/local10G/rfhorns/FlyBrain/rnaseq/analysis/v7/names_DA1_24hAPF.txt
OUTFILE_BASENAME=/local10G/rfhorns/FlyBrain/rnaseq/analysis/v7/Mz19Development_DE/df_DE_DL3_DA1_24hAPF
for i in $(seq 1 $NUM_REPLICATES)
do
    sbatch calc_DE_mannwhitneyu.sh $DF_EXPR $NAMES1 $NAMES2 $NUM_TO_SAMPLE $OUTFILE_BASENAME
done
