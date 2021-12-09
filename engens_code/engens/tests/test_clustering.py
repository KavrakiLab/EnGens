import unittest
import engens.core.FeatureSelector as fs
from engens.core.EnGens import EnGen
import pyemma
from engens.core.DimReduction import *
from engens.core.FeatureSelector import *
from engens.core.ClustEn import *

class TestDimReds(unittest.TestCase):

    def test_kmeans(self):
        test_top = "./tests/ExampleProt.pdb"
        test_traj = "./tests/ExampleTraj.xtc"
        select_expression = "residue>27 and residue<34 or residue>50 and residue<58 or residue>91 and residue<105"
        engen = EnGen(test_traj, test_top, select_expression)
        engen.init_featurizers_default()
        engen.apply_featurizations()
        feat_sele = UserFeatureSelection(2, engen)
        feat_sele.select_feature()
        reducer = dimreds["TICA"](engen)
        reducer.choose_lag(500)
        reducer.apply()
        clustering = ClusterKMeans(engen)
        params = [{"n_clusters":i} for i in range(2, 6)]
        clustering.cluster_multiple_params(params)
        clustering.plot_elbow("./tests/elbow_test.png")
        clustering.analyze_elbow_method()
        clustering.analyze_silhouette()

    def test_gmms(self):
        test_top = "./tests/ExampleProt.pdb"
        test_traj = "./tests/ExampleTraj.xtc"
        select_expression = "residue>27 and residue<34 or residue>50 and residue<58 or residue>91 and residue<105"
        engen = EnGen(test_traj, test_top, select_expression)
        engen.init_featurizers_default()
        engen.apply_featurizations()
        feat_sele = UserFeatureSelection(2, engen)
        feat_sele.select_feature()
        reducer = dimreds["TICA"](engen)
        reducer.choose_lag(500)
        reducer.apply()
        clustering = ClusterGMM(engen)
        params = [{"n_components":i} for i in range(2, 6)]
        clustering.cluster_multiple_params(params)
        clustering.plot_elbow("./tests/elbow_test.png")
        clustering.analyze_elbow_method()
        clustering.analyze_silhouette()