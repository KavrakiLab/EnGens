import unittest
import nglview as ngl
import engens.core.FeatureSelector as fs
from engens.core.EnGens import EnGen
import pyemma
import os


class TestEnGens(unittest.TestCase):


    def test_init(self):
        #correct
        test_top = "./tests/ExampleProt.pdb"
        test_traj = "./tests/ExampleTraj.xtc"
        engen = EnGen(test_traj, test_top)

        #correct with select expression
        select_expression = "residue>27 and residue<34 or residue>50 and residue<58 or residue>91 and residue<105"
        select_expression2 = [0, 1, 56, 80, 340]

        engen = EnGen(test_traj, test_top, topology_select=select_expression)
        engen = EnGen(test_traj, test_top, topology_select=select_expression2)

        #incorrect
        test_traj_fail = "./tests/jkjlk"
        test_top_fail = "./tests/kjd.xtc"
        test_traj_fail2 = 5
        test_top_fail2 = None
        select_expresion_fail = 'n238ndkjf' 

        self.assertRaises(Exception, EnGen, test_traj, test_top_fail)
        self.assertRaises(Exception, EnGen, test_traj_fail, test_top)
        self.assertRaises(Exception, EnGen, test_traj_fail, test_top_fail)
        self.assertRaises(Exception, EnGen, test_traj, test_top_fail2)
        self.assertRaises(Exception, EnGen, test_traj_fail2, test_top)
        self.assertRaises(Exception, EnGen, test_traj_fail2, test_top_fail2)
        self.assertRaises(Exception, EnGen, test_traj, test_top, select_expresion_fail)

    def test_animated_traj(self):
        test_top = "./tests/ExampleProt.pdb"
        test_traj = "./tests/ExampleTraj.xtc"
        engen = EnGen(test_traj, test_top)
        widget = engen.show_animated_traj()

        self.assertTrue(isinstance(widget, ngl.NGLWidget))

    def test_animated_traj_sele(self):
        test_top = "./tests/ExampleProt.pdb"
        test_traj = "./tests/ExampleTraj.xtc"
        select_expression = [i for i in range(50)]
        engen = EnGen(test_traj, test_top, select_expression)
        widget = engen.show_animated_traj()

        self.assertTrue(isinstance(widget, ngl.NGLWidget))

    def test_default_featurize_init(self):
        test_top = "./tests/ExampleProt.pdb"
        test_traj = "./tests/ExampleTraj.xtc"
        engen = EnGen(test_traj, test_top)
        engen.init_featurizers_default()
        self.assertEquals(len(engen.featurizers), 3)

    def test_apply_featurization(self):
        test_top = "./tests/ExampleProt.pdb"
        test_traj = "./tests/ExampleTraj.xtc"
        select_expression = [i for i in range(50)]
        default_feat2 = {
            "add_backbone_torsions": {"cossin":True, "periodic":False}
        }
        engen = EnGen(test_traj, test_top, select_expression)
        engen.add_featurizer(default_feat2)
        engen.apply_featurizations()

    def test_feat_describe(self):
        test_top = "./tests/ExampleProt.pdb"
        test_traj = "./tests/ExampleTraj.xtc"
        engen = EnGen(test_traj, test_top)
        engen.init_featurizers_default()
        res_desc = engen.describe_featurizers()
        self.assertTrue(not res_desc == "")

    def test_feat_reset(self):
        print(os.listdir("."))
        test_top = "./tests/ExampleProt.pdb"
        test_traj = "./tests/ExampleTraj.xtc"
        engen = EnGen(test_traj, test_top)
        engen.init_featurizers_default()
        engen.reset_featurizers()
        self.assertTrue(len(engen.featurizers) == 0)
        self.assertTrue(len(engen.featurizer_names) == 0)

    def test_feat_names(self):
        test_top = "./tests/ExampleProt.pdb"
        test_traj = "./tests/ExampleTraj.xtc"
        engen = EnGen(test_traj, test_top)
        engen.init_featurizers_default()
        self.assertTrue(len(engen.featurizer_names) == 3)
        print(engen.featurizer_names)
        self.assertTrue(engen.featurizer_names[0] == "residue_mindist")

    def test_choose_feat(self):
        test_top = "./tests/ExampleProt.pdb"
        test_traj = "./tests/ExampleTraj.xtc"
        engen = EnGen(test_traj, test_top)
        engen.init_featurizers_default()
        featsel = fs.UserFeatureSelection(2, engen)
        featsel.select_feature()
        featsel = fs.UserFeatureSelection(1, engen)
        featsel.select_feature()
        self.assertEquals(engen.chosen_feat_index, 1)


    def test_choose_vamp(self):
        test_top = "./tests/ExampleProt.pdb"
        test_traj = "./tests/ExampleTraj.xtc"
        select_expression = "residue>27 and residue<34 or residue>50 and residue<58 or residue>91 and residue<105"
        engen = EnGen(test_traj, test_top, select_expression)
        engen.init_featurizers_default()
        lags = [5, 10]
        dims = [2, 5]
        featsel = fs.VAMP2FeatureSelection(lags, dims, engen)
        featsel.select_feature()
        self.assertNotEqual(engen.chosen_feat_index, None)

    def test_add_pyemmafeat(self):
        test_top = "./tests/ExampleProt.pdb"
        test_traj = "./tests/ExampleTraj.xtc"
        select_expression = "residue>27 and residue<34 or residue>50 and residue<58 or residue>91 and residue<105"
        engen = EnGen(test_traj, test_top, select_expression)
        pyemma_feat = pyemma.coordinates.featurizer(engen.mdtrajref) 
        pyemma_feat.add_backbone_torsions()
        engen.add_pyemma_featurizer(pyemma_feat)
        engen.apply_featurizations()
        self.assertNotEqual(len(engen.data), 0)






        

