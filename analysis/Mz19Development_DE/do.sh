

DF_EXPR=/local10G/rfhorns/FlyBrain/rnaseq/analysis/v7/df_Mz19Development.csv
NUM_TO_SAMPLE=12
NUM_REPLICATES=100

# VA1d vs DA1

NAMES1=/local10G/rfhorns/FlyBrain/rnaseq/analysis/v7/names_VA1d_24hAPF.txt
NAMES2=/local10G/rfhorns/FlyBrain/rnaseq/analysis/v7/names_DA1_24hAPF.txt
OUTFILE_BASENAME=/local10G/rfhorns/FlyBrain/rnaseq/analysis/v7/Mz19Development_DE/df_DE_VA1d_DA1_24hAPF
for i in $(seq 1 $NUM_REPLICATES)
do
    sbatch calc_DE_mannwhitneyu.sh $DF_EXPR $NAMES1 $NAMES2 $NUM_TO_SAMPLE $OUTFILE_BASENAME
done

NAMES1=/local10G/rfhorns/FlyBrain/rnaseq/analysis/v7/names_VA1d_36hAPF.txt
NAMES2=/local10G/rfhorns/FlyBrain/rnaseq/analysis/v7/names_DA1_36hAPF.txt
OUTFILE_BASENAME=/local10G/rfhorns/FlyBrain/rnaseq/analysis/v7/Mz19Development_DE/df_DE_VA1d_DA1_36hAPF
for i in $(seq 1 $NUM_REPLICATES)
do
    sbatch calc_DE_mannwhitneyu.sh $DF_EXPR $NAMES1 $NAMES2 $NUM_TO_SAMPLE $OUTFILE_BASENAME
done

NAMES1=/local10G/rfhorns/FlyBrain/rnaseq/analysis/v7/names_VA1d_48hAPF.txt
NAMES2=/local10G/rfhorns/FlyBrain/rnaseq/analysis/v7/names_DA1_48hAPF.txt
OUTFILE_BASENAME=/local10G/rfhorns/FlyBrain/rnaseq/analysis/v7/Mz19Development_DE/df_DE_VA1d_DA1_48hAPF
for i in $(seq 1 $NUM_REPLICATES)
do
    sbatch calc_DE_mannwhitneyu.sh $DF_EXPR $NAMES1 $NAMES2 $NUM_TO_SAMPLE $OUTFILE_BASENAME
done

NAMES1=/local10G/rfhorns/FlyBrain/rnaseq/analysis/v7/names_VA1d_72hAPF.txt
NAMES2=/local10G/rfhorns/FlyBrain/rnaseq/analysis/v7/names_DA1_72hAPF.txt
OUTFILE_BASENAME=/local10G/rfhorns/FlyBrain/rnaseq/analysis/v7/Mz19Development_DE/df_DE_VA1d_DA1_72hAPF
for i in $(seq 1 $NUM_REPLICATES)
do
    sbatch calc_DE_mannwhitneyu.sh $DF_EXPR $NAMES1 $NAMES2 $NUM_TO_SAMPLE $OUTFILE_BASENAME
done

NAMES1=/local10G/rfhorns/FlyBrain/rnaseq/analysis/v7/names_adPN_adult.txt
NAMES2=/local10G/rfhorns/FlyBrain/rnaseq/analysis/v7/names_DA1_adult.txt
OUTFILE_BASENAME=/local10G/rfhorns/FlyBrain/rnaseq/analysis/v7/Mz19Development_DE/df_DE_adPN_DA1_adult
for i in $(seq 1 $NUM_REPLICATES)
do
    sbatch calc_DE_mannwhitneyu.sh $DF_EXPR $NAMES1 $NAMES2 $NUM_TO_SAMPLE $OUTFILE_BASENAME
done

# VA1d vs DC3

NAMES1=/local10G/rfhorns/FlyBrain/rnaseq/analysis/v7/names_VA1d_24hAPF.txt
NAMES2=/local10G/rfhorns/FlyBrain/rnaseq/analysis/v7/names_DC3_24hAPF.txt
OUTFILE_BASENAME=/local10G/rfhorns/FlyBrain/rnaseq/analysis/v7/Mz19Development_DE/df_DE_VA1d_DC3_24hAPF
for i in $(seq 1 $NUM_REPLICATES)
do
    sbatch calc_DE_mannwhitneyu.sh $DF_EXPR $NAMES1 $NAMES2 $NUM_TO_SAMPLE $OUTFILE_BASENAME
done

NAMES1=/local10G/rfhorns/FlyBrain/rnaseq/analysis/v7/names_VA1d_36hAPF.txt
NAMES2=/local10G/rfhorns/FlyBrain/rnaseq/analysis/v7/names_DC3_36hAPF.txt
OUTFILE_BASENAME=/local10G/rfhorns/FlyBrain/rnaseq/analysis/v7/Mz19Development_DE/df_DE_VA1d_DC3_36hAPF
for i in $(seq 1 $NUM_REPLICATES)
do
    sbatch calc_DE_mannwhitneyu.sh $DF_EXPR $NAMES1 $NAMES2 $NUM_TO_SAMPLE $OUTFILE_BASENAME
done

NAMES1=/local10G/rfhorns/FlyBrain/rnaseq/analysis/v7/names_VA1d_48hAPF.txt
NAMES2=/local10G/rfhorns/FlyBrain/rnaseq/analysis/v7/names_DC3_48hAPF.txt
OUTFILE_BASENAME=/local10G/rfhorns/FlyBrain/rnaseq/analysis/v7/Mz19Development_DE/df_DE_VA1d_DC3_48hAPF
for i in $(seq 1 $NUM_REPLICATES)
do
    sbatch calc_DE_mannwhitneyu.sh $DF_EXPR $NAMES1 $NAMES2 $NUM_TO_SAMPLE $OUTFILE_BASENAME
done

NAMES1=/local10G/rfhorns/FlyBrain/rnaseq/analysis/v7/names_VA1d_72hAPF.txt
NAMES2=/local10G/rfhorns/FlyBrain/rnaseq/analysis/v7/names_DC3_72hAPF.txt
OUTFILE_BASENAME=/local10G/rfhorns/FlyBrain/rnaseq/analysis/v7/Mz19Development_DE/df_DE_VA1d_DC3_72hAPF
for i in $(seq 1 $NUM_REPLICATES)
do
    sbatch calc_DE_mannwhitneyu.sh $DF_EXPR $NAMES1 $NAMES2 $NUM_TO_SAMPLE $OUTFILE_BASENAME
done
