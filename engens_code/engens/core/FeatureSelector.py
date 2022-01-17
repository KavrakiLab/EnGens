from engens.core.EnGens import EnGen
from typing import List
import pyemma
import warnings
import matplotlib.pyplot as plt
import numpy as np

class FeatureSelection:

    def select_feature(self) -> None:
        #implement feature selection
        pass
    

class UserFeatureSelection(FeatureSelection):

    def __init__(self, index: int, engen: EnGen) -> None:
        self.index = index
        self.engen = engen
        super().__init__()

    def select_feature(self) -> None:
        #implement feature selection
        if self.index > len(self.engen.featurizers):

            raise Exception("Featurizer index out of bounds.\
             Please chose one of the following indexes: \n" +
             self.engen.describe_featurizers())
        else:
            self.engen.chosen_feat_index = self.index  
            print("Picked featurized no. "+str(self.index)+": "+self.engen.featurizer_names[self.index])
            print(self.engen.featurizers[self.index ].describe())

class VAMP2FeatureSelection(FeatureSelection):

    def __init__(self, lags: List[int], dims: List[int], engen: EnGen) -> None:

        #dimensions and lags to try out for VAMP scoring
        self.lags = lags
        self.dims = dims
        self.scores = None
        self.engen = engen
        super().__init__()  

    def score_vamp(self, data, lag, dim):
        #do vamp scoring
        vamp = pyemma.coordinates.vamp(data, lag, dim)
        return vamp.score(data)

    def run_vamp(self) -> None:
        
        print("Choosing features with VAMP might take some time...")
        if self.engen.data == None:
            print("Generating data from featurizations")
            self.engen.apply_featurizations()
        if len(self.engen.data) == 0:
            raise Exception("No featurizers provided to EnGen.")
        
        print("Running VAMP with different parameters. Might take some time.")    
        scores = []
        for d in self.dims:
            scores_tmp = []
            for l in self.lags:
                scores_tmp_data = []
                print("dimension ={}, lag={}".format(d, l))
                for data in self.engen.data:
                    sc = self.score_vamp(data[1], l, d)
                    scores_tmp_data.append(sc.mean())
                scores_tmp.append(scores_tmp_data)
            scores.append(scores_tmp)
        self.scores = scores

    def choose_max_score_index(self) -> int:
        max_indices = {}
        for i, d in enumerate(self.dims):
            for j, l in enumerate(self.lags):
                scores_config = self.scores[i][j]
                max_val = max(scores_config)
                max_ind = scores_config.index(max_val)
                if not max_ind in max_indices:
                    max_indices[max_ind] = 1
                else:
                    max_indices[max_ind] += 1
        return max(max_indices, key=max_indices.get)

    def select_feature(self) -> None:

        print("Choosing features with VAMP might take some time...")
        if self.engen.data == None:
            print("Generating data from featurizations")
            self.engen.apply_featurizations()
        if len(self.engen.data) == 0:
            raise Exception("No featurizers provided to EnGen.")
        if len(self.engen.data) == 1:
            warnings.warn("Trying to select featurizer, only 1 provided.")
            self.engen.chosen_feat_index = 0
            return
        
        if self.scores == None:
            print("Running VAMP with different parameters.")
            self.run_vamp()
        else:
            print("Using recycled VAMP2 scores.")
        max_ind = self.choose_max_score_index()
        print("Picked featurized no. "+str(max_ind)+": "+self.engen.featurizer_names[max_ind])
        print("With maximum VAMP2 score no. "+str(max_ind))
        print(self.engen.featurizers[max_ind].describe())
        self.engen.chosen_feat_index = max_ind

    def plot_results(self) -> None:

        if self.scores == None:
            raise Exception("Can not plot results VAMP scoring has not been done yet.")

        fig, axes = plt.subplots(len(self.dims), len(self.lags), figsize=(4*len(self.lags), 3*len(self.dims)), sharey=True)
        names = None
        for i, d in enumerate(self.dims):
            for j, l in enumerate(self.lags):
                names = self.engen.featurizer_names
                for k, score in enumerate(self.scores[i][j]):
                    axes[i][j].bar(names[k], self.scores[i][j][k], label=names[k])
                axes[i][j].set_title(r'lag  $\tau$={:.1f}, dim={}'.format(l, d))
                axes[i][j].set_xticklabels(names, rotation=30)
                     
        fig.text(-0.1, 0.5, 'VAMP2 score', va='center', rotation='vertical', size=12)
        fig.tight_layout()
        handles, labels = axes[0][0].get_legend_handles_labels()
        fig.legend(handles, labels, loc='upper left')
        fig.show()


    def plot_dimensions(self, feat_ind) -> None:
        fig, ax = plt.subplots()
        for i, lag in enumerate(self.lags):
            scores = np.array([self.scores[j][i][feat_ind] for j, dim in enumerate(self.dims)])
            color = 'C{}'.format(i)
            ax.plot(self.dims, scores, '--o', color=color, label='lag={}'.format(lag))
        ax.legend()
        ax.set_xlabel('number of dimensions')
        ax.set_ylabel('VAMP2 score')
        fig.tight_layout()
        fig.show()

