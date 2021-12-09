import unittest
import nglview as ngl
from engens.core import *
from engens.core import FeatureSelector as fs
from engens.core.EnGens import EnGen
import pyemma
from engens.core.DimReduction import *
from engens.core.FeatureSelector import *

class TestDimReds(unittest.TestCase):

    def test_pca(self):
        test_top = "/home/engen/engens-code/engens/tests/ExampleProt.pdb"
        test_traj = "/home/engen/engens-code/engens/tests/ExampleTraj.xtc"
        select_expression = "residue>27 and residue<34 or residue>50 and residue<58 or residue>91 and residue<105"
        engen = EnGen(test_traj, test_top, select_expression)
        self.assertRaises(Exception, dimreds["PCA"], engen)
        engen.init_featurizers_default()
        self.assertRaises(Exception, dimreds["PCA"], engen)
        engen.apply_featurizations()
        self.assertRaises(Exception, dimreds["PCA"], engen)
        feat_sele = UserFeatureSelection(2, engen)
        feat_sele.select_feature()
        reducer = dimreds["PCA"](engen)
        reducer.plot_2d("test_pca_2d.png")
        reducer.plot_variance(90, "test_pca_var.png")
        
    def test_tica(self):
        test_top = "/home/engen/engens-code/engens/tests/ExampleProt.pdb"
        test_traj = "/home/engen/engens-code/engens/tests/ExampleTraj.xtc"
        select_expression = "residue>27 and residue<34 or residue>50 and residue<58 or residue>91 and residue<105"
        engen = EnGen(test_traj, test_top, select_expression)
        self.assertRaises(Exception, dimreds["TICA"], engen)
        engen.init_featurizers_default()
        self.assertRaises(Exception, dimreds["TICA"], engen)
        engen.apply_featurizations()
        self.assertRaises(Exception, dimreds["TICA"], engen)
        feat_sele = UserFeatureSelection(2, engen)
        feat_sele.select_feature()
        reducer = dimreds["TICA"](engen)
        reducer.plot_lag_analysis(save_loc = "test_tica_lags.png")
        reducer.choose_lag(500)
        reducer.plot_2d("test_tica_2d.png")
        reducer.plot_variance(90, "test_tica_var.png")        

    def test_hde(self):
        test_top = "/home/engen/engens-code/engens/tests/ExampleProt.pdb"
        test_traj = "/home/engen/engens-code/engens/tests/ExampleTraj.xtc"
        select_expression = "residue>27 and residue<34 or residue>50 and residue<58 or residue>91 and residue<105"
        engen = EnGen(test_traj, test_top, select_expression)
        self.assertRaises(Exception, dimreds["HDE"], engen)
        engen.init_featurizers_default()
        self.assertRaises(Exception, dimreds["HDE"], engen)
        engen.apply_featurizations()
        self.assertRaises(Exception, dimreds["HDE"], engen)
        feat_sele = UserFeatureSelection(2, engen)
        feat_sele.select_feature()
        reducer = dimreds["HDE"](engen)
        reducer.plot_lag_analysis(save_loc = "test_hde_lags.png")
        reducer.choose_lag(500)
        reducer.plot_2d("test_hde_2d.png")