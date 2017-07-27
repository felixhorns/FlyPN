from __future__ import division
import numpy as np
import scipy
import pandas as pd
import copy

import sklearn

from scipy.cluster import hierarchy

import matplotlib as mpl
from matplotlib import pyplot as plt
import seaborn as sns

import gc

def calc_DE_mannwhitneyu(X, names1, names2):
    pvalues = []
    medianA = []
    medianB = []
    meanA = []
    meanB = []
    for gene in X.index:
        A = X[names1].loc[gene]
        B = X[names2].loc[gene]
        if np.count_nonzero(A) == 0 and np.count_nonzero(B) == 0:
            pvalues.append(np.nan)
            medianA.append(0)
            medianB.append(0)
            meanA.append(0)
            meanB.append(0)
            continue
        _, pvalue = scipy.stats.mannwhitneyu(A, B)
        pvalues.append(pvalue)
        medianA.append(np.median(A))
        medianB.append(np.median(B))
        meanA.append(np.mean(A))
        meanB.append(np.mean(B))
    df_DE = pd.DataFrame({"pvalue": pvalues, "medianA": medianA, "medianB": medianB,
                          "meanA": meanA, "meanB": meanB}, index=X.index)
    df_DE.sort_values("pvalue", inplace=True)
    df_DE["pvalue_adj"] = df_DE["pvalue"] * df_DE["pvalue"].shape[0]
    return df_DE
    
def call_hits_dropout(df, dropout_cutoff=1, window_size=500, N=2000, pvalue_cutoff=0.05):
    
    df = copy.deepcopy(df)
    
    # Calculate mean
    df["mean"] = np.mean(df, axis=1)
    df.sort_values(by="mean", ascending=False, inplace=True)

    # Calculate dropouts per gene
    df["dropouts"] = np.sum(df.ix[:,:-1] < dropout_cutoff, axis=1) / df.shape[1]

    def f(df, i, window_size):
        """ Get dropout distribution for window around focal gene. 
            Calculate empirical p value. """
        
        i_upper = int(max(0, i - window_size/2))
        i_lower = int(min(df.shape[0], i + window_size/2))
        
        # Make ensemble of dropout values
        ensemble = []
        for k, (gene_name, row) in enumerate(df.iloc[i_upper:i_lower].iterrows()):
            if k + i_upper != i:
                ensemble.append(row["dropouts"])

        # Get dropouts of focal gene
        myDropouts = df.iloc[i]["dropouts"]
        
        # What fraction of ensemble is more extreme than focal gene?
        pvalue_empirical = sum(ensemble > myDropouts) / len(ensemble)
        
        return ensemble, myDropouts, pvalue_empirical

    hits = []

    for i in range(0, N):
        
        ensemble, myDropouts, pvalue_empirical = f(df, i, window_size)
        geneName = df.iloc[i].name

        if pvalue_empirical < pvalue_cutoff:
            hits.append(geneName)
    
    return hits

def call_hits_corr(df, dropout_cutoff=1, N=2000, corr_cutoff=0.8, max_canonical_vector=10, canonical_vector_spacing=3):

    df = copy.deepcopy(df)

    def f(df, C):
        """ Calculate correlations with canonical vector C """

        geneNames = []
        corrs = []
        pvalues = []

        for k, (geneName, row) in enumerate(df.iterrows()):

            myRankedExpr = row.sort_values()
            myRankedExpr = np.log10(myRankedExpr + 1)
            myY = myRankedExpr / max(myRankedExpr)
            corr, pvalue = scipy.stats.pearsonr(myY, C)
            
            geneNames.append(geneName)
            corrs.append(corr)
            pvalues.append(pvalue)

        df_corrs = pd.DataFrame({"geneName":geneNames, "corr":corrs, "pvalue":pvalues})
        
        return df_corrs

    df["mean_expr"] = np.mean(df, axis=1)
    df = df.sort_values("mean_expr", ascending=False)
    df = df.T.drop("mean_expr").T
    
    df = df.head(N)
    df[df < dropout_cutoff] = 0
    
    L = df.shape[1] # length of expr vector (num cells)

    df_corrs = pd.DataFrame()

    for k in range(1, max_canonical_vector, canonical_vector_spacing):

        C = np.zeros(L) # canonical expression profile
        C[-k:] = 1.

        df_myCorr = f(df, C)
        df_corrs = pd.concat([df_corrs, df_myCorr])

    # Filter for highly correlated hits
    df_hits = df_corrs.loc[df_corrs["corr"] > corr_cutoff].sort_values(by="corr", ascending=False)
    df_hits = df_hits.drop_duplicates(subset="geneName")
    hits = list(df_hits["geneName"])

    return hits, df_hits

def get_zscores(df, num_bin=20):

    myMean = np.mean(df, axis=1)
    myVar = np.var(df, axis=1)
    bins = np.linspace(min(myMean), max(myMean), num_bin)

    df["mean"] = myMean
    df["var"] = myVar
    df["mean_bin"] = pd.cut(myMean, bins, right=True, labels=range(1,len(bins)), include_lowest=True)

    for _, group in df.groupby("mean_bin"):

        myDispersion = np.log10(group["var"] / group["mean"])
        myDispersionStd = np.std(myDispersion)
        
        if myDispersionStd == 0: z_scores = np.zeros(len(group))

        z_scores = (myDispersion - np.mean(myDispersion)) / myDispersionStd
        df.loc[group.index, "z_score"] = z_scores

    mean = df["mean"]
    z_score = df["z_score"]
    df.drop(["mean", "var", "mean_bin", "z_score"], axis=1, inplace=True) # clean up
    
    return mean, z_score, df

def corr(df, exclude_max=0):
    # Calculate pairwise correlations between genes (columns of df)
    # Very slow because pandas corr() function is slow (compared to np.corrcoef()).
    if exclude_max > 0:
        for i in range(exclude_max):
            for col, row in enumerate(np.nanargmax(np.array(df), axis=0)):
                df.iloc[row,col] = np.nan
    return df.corr()

def get_correlated_genes(df, seeds, correlation_cutoff, min_hits=1, exclude_max=0):

    correlated_genes = []
    uncorrelated_seeds = []
    allCorrelations = pd.DataFrame(np.corrcoef(df), index=df.index, columns=df.index)

    for seed in seeds:
        
        myCorrelations = allCorrelations.loc[seed].drop(seed)
        myHits = myCorrelations[np.abs(myCorrelations) > correlation_cutoff]

        if exclude_max > 0:
            # keep only hits that are still hits after excluding maximum expressing sample            
            myValidatedHits = []
            for hit in myHits.index:
                corr_excludeMax = np.abs(corr(df.loc[[seed, hit]].T, exclude_max=exclude_max).loc[seed, hit])
                if corr_excludeMax > correlation_cutoff:
                    myValidatedHits.append(hit)
        else:
            myValidatedHits = myHits.index
        
        if len(myValidatedHits) < min_hits:
            if seed == "SpikeIn1": print "SpikeIn1 was uncorrelated seed with len(validatedHits)=", len(myValidatedHits)
            uncorrelated_seeds.append(seed)
        else:
            correlated_genes.extend(myValidatedHits)
    
    return correlated_genes, uncorrelated_seeds

def prune_singleton_correlations(df, product_cutoff=0.9):
    """ Identify genes that are driven by a single cell that highly expresses both """
    notSingletons = []
    singletons = []
    for gene1 in df.index:
        for gene2 in df.index:
            if gene1 == gene2: continue
            A = df.loc[gene1]
            B = df.loc[gene2]
            AB = A*B
            myCutoff = product_cutoff * max(A) * max(B)
            x = sum(AB > myCutoff)
            if x > 1:
                notSingletons.append(gene1)
                break
            singletons.append(gene1)
            return notSingletons, singletons

def filter_genes_overdispersed_correlates_dropouts(X, TFs, CSMs=None, exclude=None, N=50, correlation_cutoff=0.5, min_hits=1, exclude_max=0, dropout_rate_low=0.0, dropout_rate_high=1.0):
    """ Filter for informative genes for cell type identification.
        (1) Drop genes that are not expressed.
        (2) Find overdispersed genes.
        (3) Expand gene set by finding correlated genes.
        (4) Drop genes with no correlates.
        (5) Drop genes with high or low dropout rate. """
    
    # Drop genes with low max expression
    X = X.loc[np.max(X, axis=1) > 2]

    # Find overdispersed genes
    myDispersion = dispersion(X)
    myDispersion.calc_dispersion(num_bin=20) # calculate overdispersion

    hits_genome = myDispersion.get_hits(N=N)
    hits_TF = myDispersion.get_hits(N=N, candidates=TFs) # get hits among TFs
    
    if CSMs != None:
        hits_CSM = myDispersion.get_hits(N=N, candidates=CSMs) # get hits among CSMs
        hits = hits_genome.append([hits_TF, hits_CSM]).drop_duplicates().sort_values(ascending=False)
    else:
        hits = hits_genome.append([hits_TF]).drop_duplicates().sort_values(ascending=False)

    if exclude is not None:
        hits = list(set(hits.index) - set(exclude))
    else:
        hits = list(hits.index)

    # Expand gene set by finding correlated genes
    # Remove genes that have no correlates (presumably noise -- they don't belong to a "module")
    if len(hits) > 1000: print "Warning: calculating correlations between all genes and >1000 hits"
    correlated_genes, uncorrelated_seeds = get_correlated_genes(X, hits, correlation_cutoff=correlation_cutoff, min_hits=min_hits, exclude_max=exclude_max)

    hits_pruned = list(set(list(hits)) - set(list(uncorrelated_seeds)))
    hits_pruned_expanded = list(set(hits_pruned + list(set(correlated_genes))))

    Y = X.loc[hits_pruned_expanded]

    # Filter genes by dropout rate
    dropout_rate = np.sum(Y < 2, axis=1) / Y.shape[1]
    Y = Y.loc[(dropout_rate > dropout_rate_low) & (dropout_rate < dropout_rate_high)]
    
    return Y
        
class dispersion():

    def __init__(self, X):
        self.X = X
        self.X["max"] = np.max(self.X, axis=1)
        self.X.sort_values(by="max", ascending=False, inplace=True)
        self.X.drop("max", axis=1, inplace=True)
        self.mean = np.mean(self.X, axis=1)

    def calc_dispersion(self, num_bin=20):
        _, dispersion, _ = get_zscores(self.X, num_bin)
        self.dispersion = dispersion

    def plot(self, ax):
        ax.scatter(self.mean, self.dispersion)
        ax.set_xlabel("Mean expression (log2(CPM+1))")
        ax.set_ylabel("Dispersion (Z-score of log2(variance/mean))")
        ax.set_xlim(left=-0.5)

    def get_hits(self, N = None, dispersion_cutoff = None, mean_cutoff = None, candidates=None):
        if candidates != None:
            # filter genes by candidates
            dispersion_sorted = self.dispersion.loc[candidates].sort_values(ascending=False)
        else:
            dispersion_sorted = self.dispersion.sort_values(ascending=False)
        if N != None:
            # return top N hits
            hits = dispersion_sorted[:N]
        elif dispersion_cutoff != None and mean_cutoff != None:
            # return hits with dispersion and mean greater than cutoffs
            hits = dispersion_sorted.loc[self.X.loc[self.dispersion > dispersion_cutoff].loc[self.mean > mean_cutoff].index].sort_values(ascending=False)
        else:
            print "Error: N, or dispersion_cutoff and mean_cutoff must be specified"
            hits = None
        return hits
    
class PCA():

    def __init__(self, X, df, n_components):
        self.X = X
        self.df = df
        self.n_components = n_components

    def pca(self):
        self.pca = sklearn.decomposition.PCA(n_components=self.n_components)
        self.X_pca = self.pca.fit_transform(self.X.T)
        self.loadings = pd.DataFrame(self.pca.components_.T)
        self.loadings.set_index(self.X.index, inplace=True)

    def explained_variance_ratio_(self):
        return self.pca.explained_variance_ratio_

    def plot(self, ax, component_x=0, component_y=1, color_by=None):
        if color_by is not None:
            c = self.df.loc[color_by]
        else:
            c = None
        ax.scatter(self.X_pca[:,component_x], self.X_pca[:,component_y], c=c)
        ax.set_xlabel("PCA " + str(component_x + 1) + " ({0:.2%} variance)".format(self.pca.explained_variance_ratio_[component_x]))
        ax.set_ylabel("PCA " + str(component_y + 1) + " ({0:.2%} variance)".format(self.pca.explained_variance_ratio_[component_y]))
        plt.tight_layout()

    def plot_loadings(self, ax, component=0, num_genes=20):
        myLoadings = self.loadings[component].sort_values(inplace=False, ascending=False)
        plot_data = pd.concat([myLoadings.iloc[0:int(num_genes/2)], myLoadings.iloc[-int(num_genes/2):]]).sort_values(ascending=True)
        ax.barh(range(len(plot_data)), plot_data, align="center")
        ax.set_xlabel("PC " + str(component + 1) + " Loading")
        yticklabels = plot_data.index
        ax.set_yticks(range(len(plot_data)))
        ax.set_yticklabels(yticklabels)
        plt.tight_layout()

    def top_loaded_genes(self, z_score_cutoff=2, max_component=30, plot=False):
        max_component = min(self.n_components, max_component)
        hits = []
        for component in range(0, max_component):
            myLoadings = self.loadings[component].sort_values(inplace=False, ascending=False)
            myMean = np.mean(myLoadings)
            myStd = np.std(myLoadings)
            myZScores = (myLoadings - myMean) / myStd            
            myHits = list(myZScores.loc[abs(myZScores) > z_score_cutoff].index)
            hits.extend(myHits)
            if plot:
                plt.hist(myZScores, bins=np.linspace(-10, 10, 100))
                plt.xlabel("Z-score of PCA loading")
                plt.ylabel("Genes")
        hits = list(set(hits))
        return hits

class TSNE():

    def __init__(self, X, df, df_libs, n_components=2):
        self.X = X
        self.df = df
        self.n_components = n_components
        self.df_libs = df_libs

    def calc_TSNE(self, perplexity=30, random_state=0, learning_rate=500.0, early_exaggeration=4.0, metric="correlation", n_iter=1000, method="barnes_hut"):
        if metric=="correlation":
            self.dist = 1-self.X.corr()
            self.dist = np.clip(self.dist, 0.0, max(np.max(self.dist))) # clip negative values to 0 (small negative values can occur due to floating point imprecision)
            # self.D = self.dist / sum(np.sum(self.dist))
            self.D = self.dist
        elif metric=="cosine":
            self.dist = scipy.spatial.distance.squareform(scipy.spatial.distance.pdist(self.X.T, metric="cosine"))
            self.D = self.dist
        self.TSNE = sklearn.manifold.TSNE(n_components=self.n_components, metric="precomputed",
                                          perplexity=perplexity, early_exaggeration=early_exaggeration,
                                          learning_rate=learning_rate, n_iter=n_iter,
                                          method=method, verbose=2, random_state=random_state)
        self.X_tsne = self.TSNE.fit_transform(self.D)


    def plot(self, fig, ax, colorBy=None, colorMode="gene", **kwargs):
       
        if colorMode == "gene":

            if colorBy == None:
                print "Error: colorBy must be a list with 1 or 2 gene names"
                return None
            
            elif isinstance(colorBy, basestring):
                singleColor = True
                Z = self.df.loc[colorBy]
                # Z = Z / max(Z) # normalize to max
                C = Z
                
            elif len(colorBy) == 2:
                singleColor = False
                # c = [[0., 0., 0.] for _ in range(self.X_tsne.shape[0])] # initialize to black
                c = [[1., 1., 1.] for _ in range(self.X_tsne.shape[0])] # initialize to white
                color_tuple_indexes = [0, 2] # choose color channels
                for i, feature in zip(color_tuple_indexes[:len(colorBy[:2])], colorBy[:2]):
                    Z = self.df.loc[feature]
                    Z = Z / max(Z) # normalize to max
                    # Z = Z / 2 + 0.5 # compress range to 0.5 to 1 to make brighter
                    for myColor, myIntensity in zip(c, Z):
                        # myColor[i] = myIntensity
                        myColor[i] = 1. - myIntensity
                C = map(tuple, c)

            elif len(colorBy) > 2:
                print "Warning: colorBy uses a maximum of 2 colors"
                
            """
            # set uncolored points to white
            # not relevant when initialized to white
            c_whiteDefault = []
            color_intensity_cutoff = 0.1
            for color in c:
                if ((color[0] < color_intensity_cutoff) and
                    (color[1] < color_intensity_cutoff) and
                    (color[2] < color_intensity_cutoff)):
                    c_whiteDefault.append([1., 1., 1.])
                else:
                    c_whiteDefault.append(color)

            # C = map(tuple, c_whiteDefault)
            """

        elif colorMode == "genotype":
            
            C = self.df_libs.loc[self.df.columns]["color"]            

        elif colorMode == "Tissue.type":

            C = self.df_libs.loc[self.df.columns]["Tissue.type.color"]

        elif colorMode == "custom":

            C = colorBy

        if colorMode == "gene":
            norm = mpl.colors.Normalize(vmin=0, vmax=16)
            sc = ax.scatter(self.X_tsne[:,0], self.X_tsne[:,1], s=20, edgecolor="k", c=C, norm=norm, **kwargs)
        else:
            sc = ax.scatter(self.X_tsne[:,0], self.X_tsne[:,1], s=20, edgecolor="k", c=C, **kwargs)
        ax.set_xlabel("tSNE 1")
        ax.set_ylabel("tSNE 2")

        # make legend
        if colorMode == "gene" and singleColor == False:
            (left, right) = ax.get_xlim() # get current xlim
            for i, feature in zip(color_tuple_indexes, colorBy):
                c = [1., 1., 1.]
                c[i] = 0.
                ax.scatter(max(self.X_tsne[:,0]) + 1e6, 0, s=20, c=c, label=feature)
                ax.set_xlim(left, right)
                ax.legend(loc="center left", bbox_to_anchor=(1.05, 0.5)) # legend outside plot

        # make colorbar
        if colorMode == "gene" and singleColor == True:
            # cbar = fig.colorbar(sc, label="Log2(CPM+1)")
            plt.tight_layout()
            ax.set_aspect("equal")            
            # return sc, cbar
            return sc

        plt.tight_layout()
        ax.set_aspect("equal")        

        return sc

class hclust:

    def __init__(self, X, df):
        self.X = X
        self.df = df

    def cluster(self, method="average", metric="correlation"):
        pdist = scipy.spatial.distance.pdist(self.X, metric=metric)
        self.row_linkage = hierarchy.linkage(pdist, metric=metric, method=method) # cluster on gene vectors
        pdist = scipy.spatial.distance.pdist(self.X.T, metric=metric)
        pdist[np.isnan(pdist)] = 1.0
        self.col_linkage = hierarchy.linkage(pdist, metric=metric, method=method)

    def plot(self, figsize=(9,9), cmap="YlGnBu_r", **kwargs):
        cm = sns.clustermap(self.X, row_linkage=self.row_linkage, col_linkage=self.col_linkage,
                            figsize=figsize, cmap=cmap, **kwargs)
        plt.setp(cm.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
        return cm

    def get_labels(self, what=None, n_clusters=2):        
        if what == "row":
            labels = hierarchy.cut_tree(self.row_linkage, n_clusters)
        elif what == "col":
            labels = hierarchy.cut_tree(self.col_linkage, n_clusters)
        else:
            print 'Error: what must be "row" or "col"'
        return labels

    def cut(self, n_clusters_genes=2, n_clusters_cells=2):
        # Currently only works for clusters = 2
        
        gene_labels = hierarchy.cut_tree(self.row_linkage, n_clusters_genes)
        cell_labels = hierarchy.cut_tree(self.col_linkage, n_clusters_cells)

        X1 = self.X[self.X.columns[map(bool, 1-cell_labels)]]
        X2 = self.X[self.X.columns[map(bool, cell_labels)]]

        genes1 = list(self.X.T[self.X.T.columns[map(bool, 1-gene_labels)]].T.index)
        genes2 = list(self.X.T[self.X.T.columns[map(bool, gene_labels)]].T.index)

        df1 = self.df[X1.columns]
        df2 = self.df[X2.columns]

        return X1, df1, genes1, X2, df2, genes2

def get_terminal_branch_lengths(Z, max_label):
    """ Get terminal branch lengths of a linkage matrix """
    L = [d for d, label in zip(Z[:,2], Z[:,0]) if label < max_label] + [d for d, label in zip(Z[:,2], Z[:,1]) if label < max_label]
    return L

class ICIM:

    def __init__(self, X, df, TFs, CSMs, exclude, N, correlation_cutoff,
                 min_hits, exclude_max, dropout_rate_low, dropout_rate_high,
                 metric, stop_condition, N_stop, linkage_dist_stop):
        
        self.df = df
        
        self.population = {}
        self.population["0"] = X
        self.markers = {}
        self.markers["0"] = []
        self.hclust = {}
        
        # parameters
        self.TFs = TFs
        self.CSMs = CSMs
        self.exclude = exclude
        self.N = N
        self.correlation_cutoff = correlation_cutoff
        self.min_hits = min_hits
        self.exclude_max = exclude_max
        self.dropout_rate_low = dropout_rate_low
        self.dropout_rate_high = dropout_rate_high
        self.metric = metric
        self.N_stop = N_stop
        self.stop_condition = stop_condition
        self.linkage_dist_stop = linkage_dist_stop

    def step(self, parent, TFs=None, CSMs=None, exclude=None,
             N=None, correlation_cutoff=None, min_hits=None,
             exclude_max=None, dropout_rate_low=None,
             dropout_rate_high=None, metric=None, stop_condition=None,
             N_stop=None, linkage_dist_stop=None, verbose=False):
        # perform one iteration by splitting parent population

        if TFs is None: TFs = self.TFs
        if CSMs is None: CSMs = self.CSMs
        if exclude is None: exclude = self.exclude
        if N is None: N = self.N
        if correlation_cutoff is None: correlation_cutoff = self.correlation_cutoff
        if min_hits is None: min_hits = self.min_hits
        if exclude_max is None: exclude_max = self.exclude_max
        if dropout_rate_low is None: dropout_rate_low = self.dropout_rate_low
        if dropout_rate_high is None: dropout_rate_high = self.dropout_rate_high
        if metric is None: metric = self.metric
        if stop_condition is None: stop_condition = self.stop_condition
        if N_stop is None: N_stop = self.N_stop
        if linkage_dist_stop is None: linkage_dist_stop = self.linkage_dist_stop
        
        myPop = self.population[parent]
        Y = filter_genes_overdispersed_correlates_dropouts(myPop, TFs, CSMs, exclude,
                                                           N, correlation_cutoff, min_hits, exclude_max,
                                                           dropout_rate_low, dropout_rate_high)
        if verbose:
            print "Found", Y.shape[0], "genes"

        if Y.shape[0] <= 5: return []
        
        myClust = hclust(Y, self.df)
        myClust.cluster(method="average", metric=metric)
        _, X1, markerGenes1, _, X2, markerGenes2 = myClust.cut()

        self.hclust[parent] = myClust

        if stop_condition == "linkage_dist":

            Z = myClust.col_linkage
            L_terminal = get_terminal_branch_lengths(Z, myClust.X.shape[1]) # terminal branch lengths
            min_linkage_dist = min(L_terminal) # smallest terminal branch length (most similar to neighbor)

            if min_linkage_dist > linkage_dist_stop:
                # do not keep child populations and markers
                # return empty queue
                if verbose:
                    print "Failed linkage distance condition. Stopping."
                return []
        
        child1 = parent + "0"
        self.population[child1] = X1
        self.markers[child1] = markerGenes1

        child2 = parent + "1"
        self.population[child2] = X2
        self.markers[child2] = markerGenes2

        if verbose:
            print "Child populations", X1.shape[1], X2.shape[1]

        queue = []

        if stop_condition == "num_cells":
            
            if X1.shape[1] >= N_stop:
                queue.append(child1)
            
            if X2.shape[1] >= N_stop:
                queue.append(child2)
                
        elif stop_condition == "linkage_dist":

            if X1.shape[1] >= 20:
                queue.append(child1)
            if X2.shape[1] >= 20:
                queue.append(child2)

        else:

            print 'Error: stop_condition must be "num_cells" or "linkage_dist"'
            return []

        gc.collect()

        return queue

    def calc(self, skip=[], verbose=False):

        if verbose:
            print "Initial step"
        queue = []
        queue.extend(self.step("0", verbose=verbose))

        if verbose:
            print
        
        while len(queue) != 0:
            
            parent = queue.pop()
            if verbose:
                print parent
            myQueue = self.step(parent, verbose=verbose)
            myQueue = list(set(myQueue) - set(skip))
            queue.extend(myQueue)
            if verbose:
                print
            
        return None

    def get_all_markers(self):
        allMarkers = list(set([item for sublist in self.markers.values() for item in sublist]))
        return allMarkers
