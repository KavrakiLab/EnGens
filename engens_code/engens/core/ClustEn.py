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
import statistics as stat
import os
import plotly.graph_objects as go
import plotly.figure_factory as ff

import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.spatial.distance import cdist

class ClustEn(object):

    def __init__(self, engen:EnGen, clustMethod, metric, n_rep:int=5) -> None:
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
        self.thr = None
        self.n_rep = n_rep
        self.n_param_name = None
        super().__init__()

    def choose_param(self, index:int):
        self.chosen_index = index

    def choose_n(self, n:int):
        arr_n = self.get_n_params(self.params)
        p_index = np.where(arr_n==n)[0][0]
        self.choose_param(p_index)

    def get_n_params(self, params):
        params_np = np.array(params)
        return np.vectorize(lambda x: x[self.n_param_name])(params_np)

    def cluster_weights(self, i):
        pass

    def cluster_center(self, i):
        pass

    def plot_cluster_weight(self, thr=None, filename=None):
        if self.chosen_index == None: raise Exception("Choose parameters index first.")
        self.thr = thr
        weights = self.cluster_weights(self.chosen_index)
        fig = plt.figure()
        plt.bar(range(len(weights)), weights)
        plt.xlabel("Cluster number")
        plt.xticks(np.arange(0, len(weights),1))
        if not self.thr is None:
            plt.axhline(y=thr, color='r', linestyle='-')
            plt.text(0,thr+0.02, "thr={:.2f}".format(thr), color="red", ha="right", va="center")
        plt.ylabel("Cluster weight")
        if not filename == None: plt.savefig(filename)


    def plot_cluster_choice(self, filename=None):
        
        centers = self.cluster_center(self.chosen_index)
        chosen_centers = [c for i, c in enumerate(centers) if i in self.chosen_cluster_ids]
        chosen_frames = self.chosen_frames

        plt.scatter(self.data[:,0], self.data[:,1], c=self.labels[self.chosen_index],
            s=10, edgecolor='black', lw = 0.5, alpha=0.4)

        plt.scatter(np.array(chosen_centers)[:,0], np.array(chosen_centers)[:,1], c='red',
                    s=50, edgecolor='black', lw = 1)

        plt.scatter(self.data[chosen_frames][:,0], self.data[chosen_frames][:,1], c='yellow',
                    s=50, edgecolor='black', lw = 1)

        plt.legend([ "colored by cluster membership: all frames", "red: cluster_centers", "yellow: chosen_frames"],
                bbox_to_anchor=(1.05, 1.0), loc='upper left')
        plt.xlabel("C1")
        plt.ylabel("C2")
        if not filename == None: plt.savefig(filename)


    def choose_clusters(self, thr):
        if self.chosen_index == None: raise Exception("Choose parameters index first.")
        self.plot_cluster_weight(thr)
        weights = self.cluster_weights(self.chosen_index)
        clusters = [i for i, w in enumerate(weights) if w>thr]
        self.chosen_cluster_ids = clusters
        print("Chosen cluster ids: {}".format(clusters))

    def choose_conformations_center(self):
        #find cluster centers
        centers = self.cluster_center(self.chosen_index)

        closest_conf = []
        #for each chosen cluster center find the closest point within the cluster
        for c_id in self.chosen_cluster_ids:
            print(c_id)
            #extract only points within the cluster
            c_data_ind = np.argwhere(self.labels[self.chosen_index]==c_id)
            c_data = self.data[c_data_ind]
            c_data = c_data.reshape((c_data.shape[0], c_data.shape[2]))
            #find the closest point
            c_data_closest = pairwise_distances_argmin([centers[c_id]],c_data)
            #append it
            closest_conf.append(c_data_ind[c_data_closest])

        self.chosen_frames = np.array(closest_conf).flatten()

    
    def choose_conformations_hub(self):

        hub_conf = []
        #for each chosen cluster center find the closest point within the cluster
        for c_id in self.chosen_cluster_ids:
            print(c_id)
            #extract only points within the cluster
            c_data_ind = np.argwhere(self.labels[self.chosen_index]==c_id)
            c_data = self.data[c_data_ind]
            c_data = c_data.reshape((c_data.shape[0], c_data.shape[2]))
            #find the sum of the distances within the cluster
            c_data_dist = cdist(c_data,c_data)
            c_data_dist_mean = np.mean(c_data_dist.flatten())
            c_data_neigbors = np.where(c_data_dist <= c_data_dist_mean, 1, 0)
            c_data_neigbors_cnt = np.sum(c_data_neigbors, axis=0)
            c_data_hub = np.argmax(c_data_neigbors_cnt)
            #append it
            hub_conf.append(c_data_ind[c_data_hub])

        self.chosen_frames = np.array(hub_conf).flatten()

    def choose_conformations(self, filename=None, mode="center"):
        #check if parametrization is chosen 
        if self.chosen_index == None: raise Exception("Choose parameters index first.")
        
        #check if the correct mode is chosen
        if not mode in ["center", "hub"]:
            raise Exception("Mode can be: center|hub")

        #check if clusters are chosen with threshold - if not incluse all clusters in analysis
        if self.chosen_cluster_ids == None:
            n_clust = len(set(self.labels[self.chosen_index]))
            self.chosen_cluster_ids = list(range(n_clust))

        if mode == "center":
            print("Choosing points closest to the cluster center")
            self.choose_conformations_center()
        elif mode == "hub":
            print("Choosing points closest to all the other points in the cluster")
            self.choose_conformations_hub()
        self.plot_cluster_choice(filename=filename)

    def extract_conformations(self, loc:str):
        if self.chosen_frames is None: raise Exception("Choose conformations first.")

        Path(loc).mkdir(parents=True, exist_ok=True)

        for i, frameindex in enumerate(self.chosen_frames): 
            # since there was only one trajectory, trajindex will always be zero
            # if there are more trajectories, then they will be denoted by trajindex=0,1,2,...
            print("Closest conformation inside cluster " + str(self.chosen_cluster_ids[i]) + " frame " + str(frameindex) + " of the striped trajectory")
            new_file = os.path.join(str(loc),"conf_in_cluster_" + str(i)+".pdb")
            print("Extracting and saving to {}".format(new_file))
            if os.path.exists(new_file):
                print("File with the same name exists")
                print("Deleting the old file")
                os.remove(new_file)
            ret_val = call(["mdconvert -t " + self.engen.full_ref + " -o " + new_file + " -i " + str(frameindex) + " " + self.engen.full_traj_name], shell=True)
            if not ret_val == 0:
                raise(Exception("Error saving conformations."))

    def cluster_multiple_params(self, params:list):
        cls = []
        labels = []
        metric_vals = {"mean":[], "std":[]}
        self.params = params
        for p in params:
            print("Clustering with params="+str(p))
            #cluster 10 times to see the variance
            metrics = []
            for i in range(self.n_rep):
                cl = self.clustMethod(**p)
                cluster_labels = cl.fit_predict(self.data)
                metrics.append(self.compute_metric(cl))
            mean_val = metrics
            if self.n_rep > 1:
                mean_val = stat.mean(metrics)
            std_val = 0
            if self.n_rep > 1:
                std_val = stat.stdev(metrics)
                
            metric_vals["mean"].append(mean_val)
            metric_vals["std"].append(std_val)
            cls.append(cl)
            labels.append(cluster_labels)
        metric_vals["mean"] = np.array(metric_vals["mean"])
        metric_vals["std"] = np.array(metric_vals["std"])
        self.metric_vals = metric_vals
        self.labels = labels
        self.cls = cls



    # Function to plot silhouette score analysis
    def plot_silhouette(self, n_clusters, cluster_labels, X, cluster_centers, filename=None):
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
        
        if not filename == None: plt.savefig(filename)

        return silhouette_avg
    
    def cluster_center_method(self, cl, data, labels):
        pass

    def compute_metric(self, cl:GaussianMixture):
        pass

    def analyze_silhouette(self, file_prefix = None):

        avg_sil = []
        for i,p in enumerate(self.params):
            n_clusters = len(set(self.labels[i]))
            cluster_labels = self.labels[i]
            cluster_centers = self.cluster_center_method(self.cls[i], self.data, cluster_labels)
            X = self.data
            if not file_prefix is None: 
                filename = file_prefix+str(i)+".jpg"
            else:
                filename = None
            avg_sil.append(self.plot_silhouette(n_clusters, cluster_labels, X, cluster_centers, filename))
        
        best_index = avg_sil.index(max(avg_sil)) 
        print("Best parameters from silhouette analysis are "+str(self.params[best_index])) 
        return self.params[best_index]

    def plot_dendogram(self, filename=None):
        engen_dimred_distance = cdist(self.engen.dimred_data, self.engen.dimred_data, metric='euclidean')
        # get data
        data_array = engen_dimred_distance
        if self.engen.structure_names is not None:
            labels = self.engen.structure_names
        else:
            labels = self.engen.pdb_files


        # Initialize figure by creating upper dendrogram
        fig = ff.create_dendrogram(data_array, orientation='bottom', labels=labels)
        for i in range(len(fig['data'])):
            fig['data'][i]['yaxis'] = 'y2'

        # Create Side Dendrogram
        dendro_side = ff.create_dendrogram(data_array, orientation='right')
        for i in range(len(dendro_side['data'])):
            dendro_side['data'][i]['xaxis'] = 'x2'

        # Add Side Dendrogram Data to Figure
        for data in dendro_side['data']:
            fig.add_trace(data)

        # Create Heatmap
        dendro_leaves = dendro_side['layout']['yaxis']['ticktext']
        dendro_leaves = list(map(int, dendro_leaves))
        data_dist = pdist(engen_dimred_distance)
        heat_data = squareform(data_dist)
        heat_data = heat_data[dendro_leaves,:]
        heat_data = heat_data[:,dendro_leaves]

        heatmap = [
            go.Heatmap(
                x = dendro_leaves,
                y = dendro_leaves,
                z = heat_data,
                colorscale = 'Blues'
            )
        ]

        heatmap[0]['x'] = fig['layout']['xaxis']['tickvals']
        heatmap[0]['y'] = dendro_side['layout']['yaxis']['tickvals']

        # Add Heatmap Data to Figure
        for data in heatmap:
            fig.add_trace(data)

        # Edit Layout
        fig.update_layout({'width':800, 'height':800,
                                'showlegend':False, 'hovermode': 'closest',
                                })
        # Edit xaxis
        fig.update_layout(xaxis={'domain': [.15, 1],
                                        'mirror': False,
                                        'showgrid': False,
                                        'showline': False,
                                        'zeroline': False,
                                        'ticks':""})
        # Edit xaxis2
        fig.update_layout(xaxis2={'domain': [0, .15],
                                        'mirror': False,
                                        'showgrid': False,
                                        'showline': False,
                                        'zeroline': False,
                                        'showticklabels': False,
                                        'ticks':""})

        # Edit yaxis
        fig.update_layout(yaxis={'domain': [0, .85],
                                        'mirror': False,
                                        'showgrid': False,
                                        'showline': False,
                                        'zeroline': False,
                                        'showticklabels': False,
                                        'ticks': ""
                                })
        # Edit yaxis2
        fig.update_layout(yaxis2={'domain':[.825, .975],
                                        'mirror': False,
                                        'showgrid': False,
                                        'showline': False,
                                        'zeroline': False,
                                        'showticklabels': False,
                                        'ticks':""})

        # Plot!
        if filename is not None:
            fig.write_html(filename)
        return fig



class ClusterKMeans(ClustEn):

    def __init__(self, engen:EnGen, n_rep:int=5) -> None:
        super().__init__(engen, KMeans, "sum of squared distances", n_rep=n_rep)
        self.n_param_name = "n_clusters"

    def compute_metric(self, cl:KMeans):
        return cl.inertia_
    
    def cluster_center_method(self, cl, data=None, labels=None):
        return cl.cluster_centers_
    
    def cluster_weights(self, i):
        lab = self.labels[i]
        classes, counts = np.unique(lab, return_counts=True)
        return counts/sum(counts)

    def cluster_center(self, i):
        centers = self.cls[i].cluster_centers_
        return centers 


        # Function to plot elbow method analysis
    def plot_elbow(self, filename:str = None):
        x_label = [str(p["n_clusters"])for i, p in enumerate(self.params)]
        
        plt.plot(x_label, self.metric_vals["mean"], 'bx-')
        if self.n_rep > 1:
            plt.fill_between(x_label, self.metric_vals["mean"]-self.metric_vals["std"], self.metric_vals["mean"]+self.metric_vals["std"], facecolor='blue', alpha=0.4)
        plt.xlabel('number of clusters')
        plt.ylabel(self.metric)
        plt.title('Elbow Method For Optimal parameters')
        if not filename == None: plt.savefig(filename)
        plt.show()


    def analyze_elbow_method(self, filename=None):
        self.plot_elbow(filename)

        x1, y1 = 2, self.metric_vals["mean"][0]
        x2, y2 = len(self.metric_vals["mean"])+1, self.metric_vals["mean"][len(self.metric_vals["mean"])-1]

        distances = []
        for i in range(len(self.metric_vals["mean"])):
            x0 = i+2
            y0 = self.metric_vals["mean"][i]
            numerator = abs((y2-y1)*x0 - (x2-x1)*y0 + x2*y1 - y2*x1)
            denominator = sqrt((y2 - y1)**2 + (x2 - x1)**2)
            if denominator == 0:
                res = 0
            else:
                res = numerator/denominator
            distances.append(res)
        optimal_index = distances.index(max(distances))
        print("Optimal params={}".format(str(self.params[optimal_index])))
        return self.params[optimal_index]



class ClusterGMM(ClustEn):

    def __init__(self, engen:EnGen, type_ic:str="aic", n_rep:int=5) -> None:
        super().__init__(engen, GaussianMixture, type_ic, n_rep=n_rep)
        self.n_param_name = "n_components"

    def compute_metric(self, cl:GaussianMixture):
        if self.metric == "aic":
            return cl.aic(self.engen.dimred_data)
        elif self.metric == "bic":
            return cl.bic(self.engen.dimred_data)
        else:
            raise Exception("Wrong metric given to GMM!")

    def cluster_center_method(self, cl, data=None, labels=None):
        return cl.means_

    def cluster_weights(self, i):
        return self.cls[i].weights_        

    def cluster_center(self, i):
        centers = self.cls[i].means_
        return centers

    # Function to plot elbow method analysis
    def plot_ic(self, type_ic:str = "aic", filename:str = None):
        x_label = [str(p["n_components"])for i, p in enumerate(self.params)]

        plt.plot(x_label, self.metric_vals["mean"], 'bx-')
        if self.n_rep > 1:
            plt.fill_between(x_label, self.metric_vals["mean"]-self.metric_vals["std"], self.metric_vals["mean"]+self.metric_vals["std"], facecolor='blue', alpha=0.4)
        plt.xlabel('number of components')
        plt.ylabel(self.metric)
        plt.title('Information Criterion For Optimal parameters')
        if not filename == None: plt.savefig(filename)
        plt.show()

    def analyze_elbow_method(self, type_ic:str = "aic", filename:str = None):
        self.analyze_ic(type_ic, filename)
    

    def analyze_ic(self, type_ic:str = "aic", filename:str = None):

        self.plot_ic(type_ic, filename)

        x1, y1 = 2, self.metric_vals["mean"][0]
        x2, y2 = len(self.metric_vals["mean"])+1, self.metric_vals["mean"][len(self.metric_vals["mean"])-1]

        distances = []
        for i in range(len(self.metric_vals["mean"])):
            x0 = i+2
            y0 = self.metric_vals["mean"][i]
            numerator = abs((y2-y1)*x0 - (x2-x1)*y0 + x2*y1 - y2*x1)
            denominator = sqrt((y2 - y1)**2 + (x2 - x1)**2)
            distances.append(numerator/denominator)
        optimal_index = distances.index(max(distances))
        
        print("Optimal params={}".format(str(self.params[optimal_index])))

        return self.params[optimal_index]

clusterings = {
    "KM": ClusterKMeans,  
    "GMM": ClusterGMM
}