import pyemma
import sklearn
from engens.core.EnGens import EnGen
import matplotlib.pyplot as plt
import numpy as np
import plotly.express as px
from hde import HDE
import matplotlib.cm as cm
from math import sqrt
from sklearn.metrics import silhouette_samples, silhouette_score
from sklearn.cluster import KMeans
from sklearn.mixture import GaussianMixture
from sklearn.metrics import pairwise_distances_argmin
from subprocess import call
from pathlib import Path

class ClustEn(object):

    def __init__(self, engen:EnGen, clustMethod, metric) -> None:
        self.engen = engen
        self.data = engen.dimred_data
        if self.data is None: raise Exception("Data for clustering not provided")
        self.clustMethod = clustMethod
        self.metric = metric
        self.params = None
        self.cls = None
        self.labels = None
        self.metric_vals = None
        self.chosen_index = None
        self.chosen_cluster_ids = None
        self.chosen_frames = None
        super().__init__()

    def choose_param(self, index:int):
        self.chosen_index = index

    def cluster_weights(self, i):
        pass

    def cluster_center(self, i):
        pass

    def plot_cluster_weight(self, thr=None):
        if self.chosen_index == None: raise Exception("Choose parameters index first.")
        weights = self.cluster_weights(self.chosen_index)
        fig = plt.figure()
        plt.bar(range(len(weights)), weights)
        plt.xlabel("Cluster number")
        plt.xticks(np.arange(0, len(weights),1))
        if not thr is None:
            plt.axhline(y=thr, color='r', linestyle='-')
            plt.text(0,thr+0.02, "thr={:.2f}".format(thr), color="red", ha="right", va="center")
        plt.ylabel("Cluster weight")

    def plot_cluster_choice(self):
        
        centers = self.cluster_center(self.chosen_index)
        chosen_centers = [c for i, c in enumerate(centers) if i in self.chosen_cluster_ids]
        closest_conf = pairwise_distances_argmin(chosen_centers,self.data)

        plt.scatter(self.data[:,0], self.data[:,1], c=self.labels[self.chosen_index],
            s=10, edgecolor='black', lw = 0.5, alpha=0.4)

        plt.scatter(np.array(chosen_centers)[:,0], np.array(chosen_centers)[:,1], c='red',
                    s=50, edgecolor='black', lw = 1)

        plt.scatter(self.data[closest_conf][:,0], self.data[closest_conf][:,1], c='yellow',
                    s=50, edgecolor='black', lw = 1)

        plt.legend([ "colored by cluster membership: all frames", "red: cluster_centers", "yellow: chosen_frames"],
                bbox_to_anchor=(1.05, 1.0), loc='upper left')
        plt.xlabel("C1")
        plt.ylabel("C2")


    def choose_clusters(self, thr):
        if self.chosen_index == None: raise Exception("Choose parameters index first.")
        self.plot_cluster_weight(thr)
        weights = self.cluster_weights(self.chosen_index)
        clusters = [i for i, w in enumerate(weights) if w>thr]
        self.chosen_cluster_ids = clusters
        print("Chosen cluster ids: {}".format(clusters))

    def choose_conformations(self):
        if self.chosen_index == None: raise Exception("Choose parameters index first.")
        if self.chosen_cluster_ids == None:
            n_clust = len(set(self.labels[self.chosen_index]))
            self.chosen_cluster_ids = list(range(n_clust))
        centers = self.cluster_center(self.chosen_index)
        chosen_centers = [c for i, c in enumerate(centers) if i in self.chosen_cluster_ids]
        print("Chosen centers: {}".format(chosen_centers))
        closest_conf = pairwise_distances_argmin(chosen_centers,self.data)
        print("Chosen frames: {}".format(closest_conf))
        self.chosen_frames = closest_conf
        self.plot_cluster_choice()

    def extract_conformations(self, loc:str):
        if self.chosen_frames is None: raise Exception("Choose conformations first.")

        Path(loc).mkdir(parents=True, exist_ok=True)

        for i, frameindex in enumerate(self.chosen_frames): 
            # since there was only one trajectory, trajindex will always be zero
            # if there are more trajectories, then they will be denoted by trajindex=0,1,2,...
            print("Closest conformation inside cluster " + str(self.chosen_cluster_ids[i]) + " frame " + str(frameindex) + " of the striped trajectory")
            print("Extracting and saving")
            call(["mdconvert -t " + self.engen.full_ref + " -o " + str(loc) + "/conf_in_cluster_" + str(i) + ".pdb -i " + str(frameindex) + " " + self.engen.full_traj_name], shell=True)

    def cluster_multiple_params(self, params:list):
        cls = []
        labels = []
        metric_vals = []
        self.params = params
        for p in params:
            print("Clustering with params="+str(p))
            cl = self.clustMethod(**p)
            cls.append(cl)
            cluster_labels = cl.fit_predict(self.data)
            labels.append(cluster_labels)
            metric = self.compute_metric(cl)
            metric_vals.append(metric)
        self.metric_vals = metric_vals
        self.labels = labels
        self.cls = cls

    # Function to plot elbow method analysis
    def plot_elbow(self, filename:str = None):
        x_label = [str(i)for i, p in enumerate(self.params)]
        
        plt.plot(x_label, self.metric_vals, 'bx-')
        plt.xlabel('k')
        plt.ylabel(self.metric)
        plt.title('Elbow Method For Optimal parameters')
        if not filename == None: plt.savefig(filename)
        plt.show()


    def analyze_elbow_method(self):
        self.plot_elbow()

        x1, y1 = 1, self.metric_vals[0]
        x2, y2 = 100, self.metric_vals[len(self.metric_vals)-1]

        distances = []
        for i in range(len(self.metric_vals)):
            x0 = i+2
            y0 = self.metric_vals[i]
            numerator = abs((y2-y1)*x0 - (x2-x1)*y0 + x2*y1 - y2*x1)
            denominator = sqrt((y2 - y1)**2 + (x2 - x1)**2)
            distances.append(numerator/denominator)
        optimal_index = distances.index(max(distances))
        print("Optimal index={}".format(optimal_index))
        print("Optimal params={}".format(str(self.params[optimal_index])))
        return optimal_index


    # Function to plot silhouette score analysis
    def plot_silhouette(self, n_clusters, cluster_labels, X, cluster_centers):
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
    
    def cluster_center_method(self, cl, data, labels):
        pass

    def compute_metric(self, cl:GaussianMixture):
        pass

    def analyze_silhouette(self):

        avg_sil = []
        for i,p in enumerate(self.params):
            n_clusters = len(set(self.labels[i]))
            cluster_labels = self.labels[i]
            cluster_centers = self.cluster_center_method(self.cls[i], self.data, cluster_labels)
            X = self.data
            avg_sil.append(self.plot_silhouette(n_clusters, cluster_labels, X, cluster_centers))
        
        best_index = avg_sil.index(max(avg_sil)) 
        print("Best parameter index from silhouette analysis are "+str(best_index))  
        print("Best parameters from silhouette analysis are "+str(self.params[best_index])) 



class ClusterKMeans(ClustEn):

    def __init__(self, engen:EnGen) -> None:
        super().__init__(engen, KMeans, "sum of squared distances")

    def compute_metric(self, cl:GaussianMixture):
        return cl.inertia_
    
    def cluster_center_method(self, cl, data=None, labels=None):
        return cl.cluster_centers_
    
    def cluster_weights(self, i):
        lab = self.labels[i]
        classes, counts = np.unique(lab, return_counts=True)
        return counts/sum(counts)

    def cluster_center(self, i):
        return self.cls[i].cluster_centers_



class ClusterGMM(ClustEn):

    def __init__(self, engen:EnGen) -> None:
        super().__init__(engen, GaussianMixture, "bic")

    def compute_metric(self, cl:GaussianMixture):
        return cl.bic(self.engen.dimred_data)

    def cluster_center_method(self, cl, data=None, labels=None):
        return cl.means_

    def cluster_weights(self, i):
        return self.cls[i].weights_        

    def cluster_center(self, i):
        return self.cls[i].means_

clusterings = {
    "KM": ClusterKMeans,  
    "GMM": ClusterGMM
}