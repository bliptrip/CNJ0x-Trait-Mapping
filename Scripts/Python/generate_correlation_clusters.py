#!/usr/bin/env python3

import argparse
import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import cut_tree,dendrogram,inconsistent,linkage,fcluster,to_tree
import sys


# construct the argument parser and parse the arguments
def parse_args():
    parser = argparse.ArgumentParser(prog=sys.argv[0], description="Generate the distribution vectors for an annotated set of plots based on a kmeans-derived visual vocabulary of local gabor feature descriptors.")
    parser.add_argument('--correlations', type=str, default='../../Data/publication/tables/corrplot.cnj02.csv', required=False, help="CSV file containing inter-trait dependencies.")
    parsed = parser.parse_args(sys.argv[1:])
    return(parsed)


def generate_link_colors(Z, clusters):
    n = len(clusters)
    nclusters = len(np.unique(clusters))
    cmap = mpl.colormaps['tab20c'](range(nclusters+1)) #Add one for all links that are above cluster split point
    clusterids = (clusters-1).to_list() + [nclusters] * (Z.shape[0]) #One cluster
    for i in range(Z.shape[0]):
        n1 = int(Z[i,0])
        n2 = int(Z[i,1])
        if( clusterids[n1] == clusterids[n2] ):
            clusterids[i+n] = clusterids[n1]
    cmap_hex = np.apply_along_axis(lambda c: mpl.colors.to_hex(c,keep_alpha=False),1,cmap)
    return( (clusterids, cmap_hex) )

if __name__ == "__main__":
    parsed = parse_args()
    corrs = pd.read_csv(parsed.correlations, header=0, index_col=0)
    distances = 1 - np.abs(corrs)
    flatdists = squareform(distances)
    n = len(corrs)
    Z = linkage(flatdists, method='average', optimal_ordering=True) #Since our distances are measured as 1-abs(correlation), this isn't a euclidean measure and can cause issues with ward, median, 
    #R = inconsistent(Z)
    #clusters = pd.Series(fcluster(Z, 0.9), index=corrs.index, name='cluster')
    clusters = pd.Series(fcluster(Z, t=0.72, criterion='distance'), index=corrs.index, name='cluster')
    (clusterids,cmap) = generate_link_colors(Z,clusters)
    def link2color(k):
        return(cmap[clusterids[k]])
    dendrogram(Z, leaf_rotation=45, labels=corrs.index, link_color_func=link2color)
    plt.show()
    #Groupby cluster number and calculate w/in cluster inertia -- add total inertia
    print(clusters)
    
