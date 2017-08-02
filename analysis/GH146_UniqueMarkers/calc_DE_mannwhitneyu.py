from __future__ import division
import sys
import random
import numpy as np
from scipy.stats import mannwhitneyu
import pandas as pd
import uuid

def load_list(infile):
    X = []
    with open(infile) as f:
        for line in f:
            X.append(line.rstrip())
    return X

def calc_DE_mannwhitneyu(X, names1, names2):
    pvalues = []
    medianA = []
    medianB = []
    for gene in X.index:
        A = X[names1].loc[gene]
        B = X[names2].loc[gene]
        if np.count_nonzero(A) == 0 and np.count_nonzero(B) == 0:
            pvalues.append(np.nan)
            medianA.append(0)
            medianB.append(0)
            continue
        _, pvalue = mannwhitneyu(A, B)
        pvalues.append(pvalue)
        medianA.append(np.median(A))
        medianB.append(np.median(B))
    df_DE = pd.DataFrame({"pvalue": pvalues, "medianA": medianA, "medianB": medianB}, index=X.index)
    df_DE.sort_values("pvalue", inplace=True)
    df_DE["pvalue_adj"] = df_DE["pvalue"] * df_DE["pvalue"].shape[0]
    return df_DE

if __name__ == "__main__":

    infile_df_expr = sys.argv[1]
    infile_names1 = sys.argv[2]
    infile_names2 = sys.argv[3]
    outfile = sys.argv[4]

    X = pd.read_csv(infile_df_expr, header=0, index_col=0)
    names1 = load_list(infile_names1)
    names2 = load_list(infile_names2)

    df_DE = calc_DE_mannwhitneyu(X, names1, names2)

    df_DE.to_csv(outfile)

    print "Done!!"
