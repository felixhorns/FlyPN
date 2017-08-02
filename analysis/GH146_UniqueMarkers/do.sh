DF_EXPR=/local10G/rfhorns/FlyBrain/rnaseq/analysis/v7/df_GH146Alone.csv
MAX_LABEL=29

for i in $(seq 0 $MAX_LABEL)
do
    NAMES1=/local10G/rfhorns/FlyBrain/rnaseq/analysis/v7/GH146_UniqueMarkers/names_label_$i.txt
    NAMES2=/local10G/rfhorns/FlyBrain/rnaseq/analysis/v7/GH146_UniqueMarkers/names_notLabel_$i.txt
    OUTFILE_BASENAME=/local10G/rfhorns/FlyBrain/rnaseq/analysis/v7/GH146_UniqueMarkers/df_DE_GH146_UniqueMarkers_label_$i.csv
    sbatch calc_DE_mannwhitneyu.sh $DF_EXPR $NAMES1 $NAMES2 $OUTFILE_BASENAME
done
