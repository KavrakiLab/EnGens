import unittest
import nglview as ngl
from engens.core.EnGens import EnGen
import pyemma
from engens.core.DimReduction import *
from engens.core.FeatureSelector import *
from engens.core.ClustEn import *
from engens.core.PlotUtils import PlotUtils
import os



class TestPlots(unittest.TestCase):

    def __init__(self, *args, **kwargs) -> None:
        
        super(TestPlots, self).__init__(*args, **kwargs)
        test_top = "./tests/ExampleProt.pdb"
        test_traj = "./tests/ExampleTraj.xtc"
        select_expression = "residue>1 and residue<9 or residue>50 and residue<58 or residue>91 and residue<105"
        engen = EnGen(test_traj, test_top, select_expression)
        engen.init_featurizers_default()
        engen.apply_featurizations()
        feat_sele = UserFeatureSelection(1, engen)
        feat_sele.select_feature()
        reducer = dimreds["TICA"](engen)
        reducer.choose_lag(10)
        reducer.apply()
        self.engen = engen
        self.feat_sele = feat_sele  
        self.dimred = reducer      
        clustering = ClusterKMeans(self.engen, n_rep=3)
        params = [{"n_clusters":i} for i in range(6,8)]
        clustering.cluster_multiple_params(params)
        clustering.choose_param(1)
        clustering.choose_conformations()
        clustering.extract_conformations(".")
        self.clust = clustering

    def test_init(self):
        p_utils1 = PlotUtils(self.engen, self.clust)
        p_utils2 = PlotUtils(self.engen, self.clust, "./tests/plot_util_plots")

    def test_db1(self):
        p_utils1 = PlotUtils(self.engen, self.clust)
        p_utils1.dashboard1()
        p_utils2 = PlotUtils(self.engen, self.clust, "./tests/plot_util_plots")
        p_utils2.dashboard1()

