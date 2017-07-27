from __future__ import division
import sys
import random
import numpy as np
import scipy
import pandas as pd
import uuid

def load_list(infile):
    X = []
    with open(infile) as f:
        for line in f:
            X.append(line.rstrip())
    return X

# Define functions for information theory analysis

# Compute mutual information between each gene and cluster assignment
from sklearn.metrics import mutual_info_score
def calc_mutual_information(df_discrete, labels):
    mutual_informations = []
    for symbol, row in df_discrete.iterrows():
        mi = mutual_info_score(row, labels) / np.log(2) # calculate mutual information, convert from nats to bits
        mutual_informations.append(mi)
    df_result = pd.DataFrame()
    df_result["symbol"] = df_discrete.index
    df_result.set_index("symbol", inplace=True)
    df_result["mutual_information"] = mutual_informations
    df_result.sort_values("mutual_information", inplace=True, ascending=False)
    return df_result

# Calculate mean expression of each gene within each cluster
def calc_summary_by_label(df, labels, summary_func=np.mean):
    df_temp = copy.deepcopy(df.T)
    df_temp["label"] = labels
    df_summary_by_label = df_temp.groupby("label").aggregate(summary_func).T
    del df_temp
    return df_summary_by_label

# Calculate total information of a gene set
def calc_cumulative_information_of_set(df_discrete, genes, df_labels):
    """ Calculates total information of gene set """
    # Entropy of class without knowledge of genes
    H_naive = scipy.stats.entropy(df_labels["label"].value_counts(normalize=True), base=2)
    # Get discretized gene expression matrix with only selected genes
    df_temp = df_discrete.loc[genes].T

    H = 0
    # Find unique combinations of expression levels
    for _, row in df_temp.drop_duplicates().iterrows():
        cell_names = df_temp.index[np.all(df_temp == row, axis=1)] # Get cell names having this unique combination of expression levels
        labels_cond = df_labels.loc[cell_names]["label"] # Get class labels of cells
        # Calculate entropy of classification (conditional on expression levels)
        H_cond = scipy.stats.entropy(labels_cond.value_counts(normalize=True), base=2)
        weight = len(cell_names) / df_temp.shape[0] # Weight by fraction of cells in this set
        H += H_cond*weight

    I = H_naive - H
    return I

# Calculate total information for each gene set defined by iteratively adding a single gene
# starting from the beginning of the list genes
def calc_cumulative_informations(df_discrete, genes, df_labels, N=5):
    """ Calculates total information of gene sets starting from top of df """
    cis = []
    for i in range(0,N):
        my_ci = calc_cumulative_information_of_set(df_discrete, genes[:i+1], df_labels)
        cis.append(my_ci)
    cis = np.array(cis)
    df_result = pd.DataFrame()
    df_result["symbol"] = genes
    df_result.set_index("symbol", inplace=True)
    df_result["cumulative_mutual_information"] = ""
    df_result["cumulative_mutual_information"] = np.nan
    df_result["cumulative_mutual_information"].loc[genes[:N]] = cis
    return df_result

# Find nonredundant gene set
# Define function to pick non-redundant genes

def find_nonredundant_gene_set(df_discrete, genes,
                               df_labels,
                               df_info,
                               H_naive,
                               N_constrain=20,
                               cumulative_information_cutoff=0.99,
                               verbose=False):
        
    cis = []
    nonredundant_genes = []
    current_ci = 0
    current_relative_ci = 0.0
    
    # Rank genes by mutual information
    remaining_genes = list(df_info.loc[genes]["mutual_information"].sort_values(ascending=False).head(n=N_constrain).index)
    
    i = 0
    
    while True:
        
        i +=1
        
        # Calculate information gain for each gene
        df_info_gains = pd.DataFrame(index=list(remaining_genes))
        my_info_gains = []
        
        for gene in remaining_genes:
            my_genes = nonredundant_genes + [gene]
            my_ci = calc_cumulative_information_of_set(df_discrete, my_genes, df_labels)
            my_info_gain = my_ci - current_ci
            my_info_gains.append(my_info_gain)
        
        df_info_gains["info_gain"] = my_info_gains
        
        # Sort genes by information gain
        df_info_gains.sort_values("info_gain", ascending=False, inplace=True)
        
        # Take best gene
        hit = df_info_gains["info_gain"].index[0]
        
        nonredundant_genes.append(hit)
        remaining_genes.remove(hit)
        current_ci = current_ci + df_info_gains.iloc[0]["info_gain"]
        current_relative_ci = current_ci / H_naive
        
        if verbose:
            print hit
            print current_relative_ci
            print
        
        if len(remaining_genes) == 0 or current_relative_ci > 0.99:
            return nonredundant_genes

# Load data
infile_genes = "/local10G/rfhorns/FlyBrain/rnaseq/analysis/v7/genes_genome_noTFs_noCSMs_v2.txt"
infile_df_info = "/local10G/rfhorns/FlyBrain/rnaseq/analysis/v7/df_info.csv"
infile_df_discrete = "/local10G/rfhorns/FlyBrain/rnaseq/analysis/v7/df_discrete.csv"
infile_df_labels = "/local10G/rfhorns/FlyBrain/rnaseq/analysis/v7/labels_HDBSCAN.csv"

genes_genome_noTFs_noCSMs = load_list(infile_genes)
df_info = pd.read_csv(infile_df_info, index_col=0, header=0)
df_discrete = pd.read_csv(infile_df_discrete, index_col=0, header=0)
df_labels = pd.read_csv(infile_df_labels, sep="\t", index_col=0, header=0)

H_naive = 4.51456791675 # calculated in notebook

# Define parameters
N_genes = 1000
outfile_basename = "df_info_nonredundant_random"

# Select random genes
random_genes = np.random.choice(genes_genome_noTFs_noCSMs, N_genes, replace=False)
random_genes_topInfo = df_info.loc[random_genes].sort_values("mutual_information", ascending=False).head(n=30).index

genes_nonredundant_randomGenes = find_nonredundant_gene_set(df_discrete, random_genes_topInfo,
                                                            df_labels, df_info, H_naive,
                                                            N_constrain=30,
                                                            cumulative_information_cutoff=0.99,
                                                            verbose=False)

df_info_nonredundant_random = df_info.copy()

# Calculate cumulative information for top N genes
cumulative_informations = calc_cumulative_informations(df_discrete, genes_nonredundant_randomGenes, df_labels, N=len(genes_nonredundant_randomGenes))
df_info_nonredundant_random["cumulative_information"] = cumulative_informations["cumulative_mutual_information"]

# Calculate information relative to total entropy
relative_cumulative_informations = df_info_nonredundant_random["cumulative_information"] / H_naive
df_info_nonredundant_random["relative_cumulative_information"] = relative_cumulative_informations

# Write output df to file
suffix = str(uuid.uuid4())
outfile = outfile_basename + "." + suffix + ".csv"
df_info_nonredundant_random.loc[genes_nonredundant_randomGenes].to_csv(outfile)

print "Done!!"
