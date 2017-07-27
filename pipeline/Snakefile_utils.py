import os
import re
import subprocess

##### Function for loading seeds
def load_seeds(infile):
    seeds = []
    with open(infile, 'rU') as f:
        for line in f:
            seeds.append(line.rstrip())
    return seeds

##### Functions for transferring files to/from cluster
def name_on_scratch(s, scratch):
    return scratch+"/"+os.path.basename(s)

def names_on_scratch(names, scratch):
    return [name_on_scratch(n, scratch) for n in names]

def cp_to_scratch(inputs, scratch):
    for i in inputs:
      cmd = "rsync -aW " + i + " " + name_on_scratch(i, scratch)
      subprocess.call(cmd, shell=True)
    return None

def cp_from_scratch(outputs, scratch):
    for o in outputs:
    	cmd = "rsync -aW " + name_on_scratch(o, scratch) + " " + o
	subprocess.call(cmd, shell=True)
    return None

##### Functions for getting file names and unzipping files

def get_all_files(d):
    return [d+"/"+f for f in os.listdir(d) if os.path.isfile(os.path.join(d, f))]

def unzip_fastq(f):
    cmd = "gunzip " + f
    p = subprocess.Popen(cmd, shell=True)
    return None

def get_raw_fastqs(wildcards):

    # unzip all fastqs
    all_files = get_all_files(wildcards.dir)
    fastq_gzs = [f for f in all_files if ".fastq.gz" in f]
    for f in fastq_gzs:	unzip_fastq(f)

    # find the raw fastqs
    all_files = get_all_files(wildcards.dir)
    raw_fastqs = []
    for f in all_files:
        if ("L00" in f) and ("R"+wildcards.R in f) and (os.path.splitext(f)[1] == ".fastq"):
            raw_fastqs.append(f)
    return raw_fastqs

def get_all_fastq_gzs_R1(wildcards):
    all_files = get_all_files(wildcards.dir)
    fastq_gzs = []
    for f in all_files:
        if ("_R1_" in f) and (".fastq" in f) and (os.path.splitext(f)[1] == ".gz"):
            fastq_gzs.append(f)
    return sorted(fastq_gzs)

def get_all_fastq_gzs_R2(wildcards):
    all_files = get_all_files(wildcards.dir)
    fastq_gzs = []
    for f in all_files:
        if ("_R2_" in f) and (".fastq" in f) and (os.path.splitext(f)[1] == ".gz"):
            fastq_gzs.append(f)
    return sorted(fastq_gzs)
