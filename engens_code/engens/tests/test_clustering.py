import unittest
import engens.core.FeatureSelector as fs
from engens.core.EnGens import EnGen
import pyemma
from engens.core.DimReduction import *
from engens.core.FeatureSelector import *
from engens.core.ClustEn import *
import os

class TestClust(unittest.TestCase):

    def __init__(self, *args, **kwargs) -> None:
        
        super(TestClust, self).__init__(*args, **kwargs)
        test_top = "./tests/ExampleProt.pdb"
        test_traj = "./tests/ExampleTraj.xtc"
        select_expression = "residue>1 and residue<9 or residue>50 and residue<58 or residue>91 and residue<105"
        engen = EnGen(test_traj, test_top, select_expression)
        engen.init_featurizers_default()
        engen.apply_featurizations()
        feat_sele = UserFeatureSelection(2, engen)
        feat_sele.select_feature()
        reducer = dimreds["TICA"](engen)
        reducer.choose_lag(500)
        reducer.apply()
        self.engen = engen

    def test_kmeans(self):
        clustering = ClusterKMeans(self.engen, n_rep=3)
        params = [{"n_clusters":i} for i in range(2, 6)]
        clustering.cluster_multiple_params(params)
        clustering.plot_elbow("./tests/elbow_test.png")
        clustering.analyze_elbow_method()
        clustering.analyze_silhouette()

    def test_extracting_closest_confs_KM(self):
        #test KMeans
        clustering = ClusterKMeans(self.engen, n_rep=3)
        params = [{"n_clusters":i} for i in range(2, 6)]
        clustering.cluster_multiple_params(params)
        clustering.choose_param(3)
        clustering.choose_conformations()
        clustering.extract_conformations(".")
        for i, elem in enumerate(clustering.chosen_frames):
            self.assertTrue(clustering.labels[clustering.chosen_param_index][elem] == clustering.chosen_cluster_ids[i])
        #test GMM
        clustering = ClusterGMM(self.engen, n_rep=3)
        params = [{"n_components":i} for i in range(2, 6)]
        clustering.cluster_multiple_params(params)
        clustering.choose_param(3)
        clustering.choose_conformations()
        clustering.extract_conformations(".")
        for i, elem in enumerate(clustering.chosen_frames):
            self.assertTrue(clustering.labels[clustering.chosen_param_index][elem] == clustering.chosen_cluster_ids[i])
        


    def test_gmms(self):

        clustering = ClusterGMM(self.engen,n_rep=3)
        params = [{"n_components":i} for i in range(2, 6)]
        clustering.cluster_multiple_params(params)
        clustering.plot_ic(filename = "./tests/aic_test.png")
        clustering.analyze_ic()
        clustering = ClusterGMM(self.engen, type_ic="bic")
        params = [{"n_components":i} for i in range(2, 6)]
        clustering.cluster_multiple_params(params)
        clustering.plot_ic(filename = "./tests/bic_test.png")
        clustering.analyze_ic()
        clustering.analyze_silhouette()