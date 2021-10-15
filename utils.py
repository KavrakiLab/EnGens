import mdtraj as md
import numpy as np
import pyemma
import matplotlib.pyplot as plt
import itertools
from sklearn.cluster import KMeans
from scipy.spatial.distance import cdist, pdist
from math import sqrt
import wget
import sys
import os
from datetime import date
from subprocess import call
from mpl_toolkits.mplot3d import Axes3D 
import glob
from pathlib import Path
from Bio.PDB.PDBParser import PDBParser
from sklearn.datasets import make_blobs
from sklearn.metrics import silhouette_samples, silhouette_score
import matplotlib.cm as cm
import json

#util function to identify indexes of given residues
def extract_res_ind(top_file, res_df):
    # Use biopython to identify the chains and their size
    parser = PDBParser(PERMISSIVE=1)
    structure = parser.get_structure("struct", top_file)
    residues = [elem for elem in structure.get_residues()]
    chains = [c for c in structure.get_chains()]
    chain_start_index = [0]+[ len(chains[i-1]) for i in range(1, len(chains))]
    chain_res = {}
    chain_sind = dict(zip([c.id for c in chains], chain_start_index))
    for c in chains:
        res = [elem for elem in c.get_residues()]
        chain_res[c.id] = res
    results = []
    for ind, row in res_df.iterrows():
        c = row["chain"]
        res_id = row["residue_id"]
        for i, elem in enumerate(chain_res[c]):
            if res_id==elem.get_resname()+str(elem.get_full_id()[3][1]):
                index = chain_sind[c]+i
                results.append(index)
    return results
    

# util function to find a residue segment in a biopython structure 
def lookup_chain_cdr(chain, cdr):
    for i, residue in enumerate(chain[:-len(cdr)]):
        part = chain[i:i+len(cdr)]
        part_aa = [elem.get_resname() for elem in part]
        if part_aa == cdr:
            return ([elem.get_resname()+str(elem.get_full_id()[3][1]) for elem in part], [j for j in range(i, i+len(cdr))])
    return ([], None)

#Function to calculate wcss
def calculate_wcss(data):
    wcss = []
    for n in range(1, 100):
        kmeans = KMeans(n_clusters=n)
        kmeans.fit(X=data)
        wcss.append(kmeans.inertia_)
    return wcss
    
# Function to define the optimal number of clusters
def optimal_number_of_clusters(wcss):
    x1, y1 = 1, wcss[0]
    x2, y2 = 100, wcss[len(wcss)-1]

    distances = []
    for i in range(len(wcss)):
        x0 = i+2
        y0 = wcss[i]
        numerator = abs((y2-y1)*x0 - (x2-x1)*y0 + x2*y1 - y2*x1)
        denominator = sqrt((y2 - y1)**2 + (x2 - x1)**2)
        distances.append(numerator/denominator)
    
    return distances.index(max(distances)) + 2

# Function to plot elbow method analysis
def plot_elbow(K, ssd):
    plt.plot(K, ssd, 'bx-')
    plt.xlabel('k')
    plt.ylabel('Sum_of_squared_distances')
    plt.title('Elbow Method For Optimal k')
    plt.savefig("kmeans_elbow.png")
    plt.show()

# Function to plot silhouette score analysis
def plot_silhouette(n_clusters, cluster_labels, X, cluster_centers):
    fig, (ax1, ax2) = plt.subplots(1, 2)
    fig.set_size_inches(18, 7)

    # The 1st subplot is the silhouette plot
    # The silhouette coefficient can range from -1, 1 but in this example all
    # lie within [-0.1, 1]
    ax1.set_xlim([-0.1, 1])
    # The (n_clusters+1)*10 is for inserting blank space between silhouette
    # plots of individual clusters, to demarcate them clearly.
    ax1.set_ylim([0, len(X) + (n_clusters + 1) * 10])


    # The silhouette_score gives the average value for all the samples.
    # This gives a perspective into the density and separation of the formed
    # clusters
    silhouette_avg = silhouette_score(X, cluster_labels)
    
    # Compute the silhouette scores for each sample
    sample_silhouette_values = silhouette_samples(X, cluster_labels)

    y_lower = 10
    for i in range(n_clusters):
        # Aggregate the silhouette scores for samples belonging to
        # cluster i, and sort them
        ith_cluster_silhouette_values = \
            sample_silhouette_values[cluster_labels == i]

        ith_cluster_silhouette_values.sort()

        size_cluster_i = ith_cluster_silhouette_values.shape[0]
        y_upper = y_lower + size_cluster_i

        color = cm.nipy_spectral(float(i) / n_clusters)
        ax1.fill_betweenx(np.arange(y_lower, y_upper),
                          0, ith_cluster_silhouette_values,
                          facecolor=color, edgecolor=color, alpha=0.7)

        # Label the silhouette plots with their cluster numbers at the middle
        ax1.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i))

        # Compute the new y_lower for next plot
        y_lower = y_upper + 10  # 10 for the 0 samples

    ax1.set_title("The silhouette plot for n_clusters="+str(n_clusters)+" avg score is "+str(silhouette_avg))
    ax1.set_xlabel("The silhouette coefficient values")
    ax1.set_ylabel("Cluster label")

    # The vertical line for average silhouette score of all the values
    ax1.axvline(x=silhouette_avg, color="red", linestyle="--")

    ax1.set_yticks([])  # Clear the yaxis labels / ticks
    ax1.set_xticks([-0.1, 0, 0.2, 0.4, 0.6, 0.8, 1])

    # 2nd Plot showing the actual clusters formed
    colors = cm.nipy_spectral(cluster_labels.astype(float) / n_clusters)
    ax2.scatter(X[:, 0], X[:, 1], marker='.', s=30, lw=0, alpha=0.7,
                c=colors, edgecolor='k')

    # Labeling the clusters
    centers = cluster_centers
    # Draw white circles at cluster centers
    ax2.scatter(centers[:, 0], centers[:, 1], marker='o',
                c="white", alpha=1, s=200, edgecolor='k')

    for i, c in enumerate(centers):
        ax2.scatter(c[0], c[1], marker='$%d$' % i, alpha=1,
                    s=50, edgecolor='k')

    ax2.set_title("The visualization of the clustered data.")
    ax2.set_xlabel("Feature space for the 1st feature")
    ax2.set_ylabel("Feature space for the 2nd feature")

    plt.suptitle(("Silhouette analysis for KMeans clustering on sample data "
                  "with n_clusters = %d" % n_clusters),
                 fontsize=14, fontweight='bold')
    return silhouette_avg

# make restraints for piper
#inputs used for restraints_MMR.json
#peptide_res=['C1', 'C4', 'C5', 'C6', 'C6', 'C6']
#cdr_res=["A98", "B98", "B99", "B99", "B97", "B96"]
#range_min=[2, 2, 3, 3, 2, 4]
#range_max=[8, 8, 9, 9, 8, 10]
def make_restraints(peptide_res, cdr_res, range_min, range_max, res_fname):
    restraints = {}
    restraints["groups"] = []
    restraints["groups"].append({})
    restraints["groups"][0]["restraints"] = []
    restraints["groups"][0]["required"] = 1
    for pep_elem, cdr_elem, rmin, rmax in zip(peptide_res, cdr_res, range_min, range_max):
        res = {}
        res["rec_resid"]=pep_elem[1:]
        res["rec_chain"]=pep_elem[0]
        res["lig_resid"]=cdr_elem[1:]
        res["lig_chain"]=cdr_elem[0]
        res["dmax"]=rmax
        res["dmin"]=rmin
        res["type"]="residue"
        restraints["groups"][0]["restraints"].append(res)
    with open(res_fname, 'w') as fp:
        json.dump(restraints, fp)