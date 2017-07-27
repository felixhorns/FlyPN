#!/bin/bash

source /local10G/rfhorns/resources/anaconda2/envs/py35/bin/activate /local10G/rfhorns/resources/anaconda2/envs/py35/

MY_HOME=/local10G/rfhorns/FlyBrain/rnaseq
CONFIGFILE=$MY_HOME/scripts/config.json

SNAKEFILE=$MY_HOME/scripts/Snakefile.py

# Specify log file
DATETIME=$(date "+%Y_%m_%d_%H_%M_%S")
LOGFILE=$MY_HOME/log/Snakefile.$DATETIME.log

# Unlock dir
# snakemake all --snakefile $SNAKEFILE --configfile $CONFIGFILE --cluster "sbatch --ntasks=1 --job-name={params.name} --cpus-per-task={threads} --partition={params.partition} --mem={params.mem} -o {params.name}.%j.log" --keep-target-files -j 200 -w 100 -k -r -n --rerun-incomplete --unlock

# Dry run snakemake
# snakemake all --snakefile $SNAKEFILE --configfile $CONFIGFILE --cluster "sbatch --ntasks=1 --job-name={params.name} --cpus-per-task={threads} --partition={params.partition} --mem={params.mem} -o {params.name}.%j.log" --keep-target-files -j 200 -w 100 -k -r -n --rerun-incomplete

# Run snakemake
nohup snakemake all --snakefile $SNAKEFILE --configfile $CONFIGFILE --cluster "sbatch --ntasks=1 --job-name={params.name} --cpus-per-task={threads} --partition={params.partition} --mem={params.mem} -o {params.name}.%j.log" --keep-target-files -j 200 -w 100 -k --rerun-incomplete > $LOGFILE &

echo Log is
echo $LOGFILE
echo
