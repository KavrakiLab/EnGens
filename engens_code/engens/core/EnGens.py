
import pytraj as pt
import nglview as ngl
import mdtraj
import pyemma.coordinates 
import pyemma.coordinates.data
import os

'''
EnGen - the main class for ensemble generation and analysis 
'''

class EnGen(object):

    def __init__(self, trajectory, topology, topology_select=None):
        """ generates and visualizes conformational ensembles for ensemble docking
        
        Parameters
        ------------
        trajectory: str 
            a path to the trajectory file (.xtc, .dcd, etc.)
        topology: str
            a path to the trajectory file (.pdb)
        topology_select: str | list, default = None
            a select string that indicates what part of the topology is relevant for featurization
            or an array of atom indices of interest (as required by mdtraj)
        """

        self.traj = trajectory
        self.traj_name = trajectory
        self.full_traj_name = trajectory
        self.ref = topology
        self.full_ref = topology
        self.refrestraint = topology_select
        self._selection_indices = None

        # load a restrained trajectory based on a selection/list/None
        self.mdtrajref = mdtraj.load(self.ref).topology
        if isinstance(self.refrestraint, str):
            tmp_top = mdtraj.load(self.ref).topology
            self._selection_indices = tmp_top.select(topology_select)
            self.mdtrajref = mdtraj.load(self.ref, atom_indices=self._selection_indices).topology
            tmp_traj = mdtraj.load(trajectory, top=topology, atom_indices=self._selection_indices)
            tmp_name = trajectory[:trajectory.rfind(".")]+"-engen-selected.xtc"
            tmp_pdb_name = trajectory[:trajectory.rfind(".")]+"-engen-selected.pdb"
            tmp_traj.save(tmp_name)
            tmp_traj[0].save(tmp_pdb_name)
            self.traj = tmp_name
            self.traj_name = tmp_name
            self.ref = tmp_pdb_name

        elif isinstance(self.refrestraint, list):
            self._selection_indices = self.refrestraint
            self.mdtrajref = mdtraj.load(self.ref, atom_indices=self._selection_indices).topology
            tmp_traj = mdtraj.load(trajectory, top=topology, atom_indices=self._selection_indices)
            tmp_name = trajectory[:trajectory.rfind(".")]+"-engen-selected.xtc"
            tmp_pdb_name = trajectory[:trajectory.rfind(".")]+"-engen-selected.pdb"
            tmp_traj.save(tmp_name)
            tmp_traj[0].save(tmp_pdb_name)
            self.traj = tmp_name
            self.traj_name = tmp_name
            self.ref = tmp_pdb_name


        self.featurizers = []
        self.featurizer_names =[]
        self.data = None
        self.dimRed = None
        self.cluster = None 
        self.chosen_feat_index = -1
        self.dimred_data = None

    @property
    def traj_name(self):
        return self._traj_name

    @traj_name.setter
    def traj_name(self, t):
        valid_ext =[".xtc", ".dcd", ".trr", ".binpos", ".netcdf", ".arc", ".hdf5", ".lammpstrj"]
        if not t: raise Exception("Trajectory can not be empty")
        if not (isinstance(t, str) or isinstance(t, list)): raise Exception("Trajectory must be a string")
        if not t[t.rfind('.'):] in valid_ext: raise Exception("Trajectory must be with a correct extension")
        if not os.path.exists(t): Exception("Trajectory does not exist.")
        self._traj_name = t 

    @property
    def ref(self):
        return self._ref

    @ref.setter
    def ref(self, r):
        valid_ext =[".pdb"]
        if not r: raise Exception("Topology reference can not be empty")
        if not isinstance(r, str): raise Exception("Topology reference must be a string")
        if not r[r.rfind('.'):] in valid_ext: raise Exception("Topology reference must be with a correct extension")
        if not os.path.exists(r): Exception("Topology reference does not exist.")
        self._ref = r
    
    @property
    def refrestraint(self):
        return self._refrestraint

    @refrestraint.setter
    def refrestraint(self, r):
        if not r: 
            self._refrestraint = r
            return
        if not (isinstance(r, str) or isinstance(r, list)): 
            raise Exception("Topology restraint must be a string")
        self._refrestraint = r
    
    def show_animated_traj(self):
        """
        Returns an animated nglview widget with the loaded trajectory

        Returns
        -----------
        widget: NGLWidget object with loaded pytraj

        """
        pt_traj = pt.load(self.traj_name, self.ref)
        widget = ngl.show_pytraj(pt_traj, gui=True)
        return widget

    #--------------------FEATURIZATION--------------------------#
    
    def reset_featurizers(self):
        self.featurizer_names = []
        self.featurizers = []

    def add_featurizer(self, feats: dict):
        """
        Adds another featurization type to the list of featurizers
        
        Parameters
        ------------
        feats: dict object with entries {"add_feature_func_name": params, ...}
                "add_feature_func_name" should be an add function of pyemma.coordinates.featurizer.MDFeaturizer
                params should be parameters of the given function
        """
        pyemma_feat = pyemma.coordinates.featurizer(self.mdtrajref) 
        name = ""
        for key, params in feats.items():
            name+=key[len("add_"):]
            name+="&"
            func =  getattr(pyemma_feat, key)
            func(**params)
        name = name[:-1]
        self.featurizer_names.append(name)
        self.featurizers.append(pyemma_feat)

    def add_pyemma_featurizer(self, pyemma_feat: pyemma.coordinates.data.MDFeaturizer, name: str):
        """
        Adds another featurization type to the list of featurizers
        
        Parameters
        ------------
        feat: pyEmma featurizer
        name: name for the featurizer
        """
        self.featurizer_names.append(name)
        self.featurizers.append(pyemma_feat)

    def init_featurizers_default(self):
        """
        Adds default featurizers type to the list of featurizers
        
        Parameters
        ------------
        type: type of function for adding to inital default featurizer

        """

        # only called if initialization
        if not len(self.featurizers) == 0:
            return

        default_feat1 = {
            "add_residue_mindist": {"scheme":'closest-heavy'}
        }
        default_feat2 = {
            "add_backbone_torsions": {"cossin":True, "periodic":False}

        }
        default_feat3 = {
            "add_backbone_torsions": {"cossin":True, "periodic":False},
            "add_residue_mindist": {"scheme":'closest-heavy'}

        }
        default_feats = [default_feat1, default_feat2, default_feat3]

        for feat_dict in default_feats:
            self.add_featurizer(feat_dict)
        
    def apply_featurizations(self):
        self.data = []
        for f in self.featurizers:
            self.data+=pyemma.coordinates.source(self.traj, features=f)

    def describe_featurizers(self):
        res = ""
        for i, f in enumerate(self.featurizers):
            res += "Featurizer no. "+str(i)+":\n "
            res += self.featurizer_names[i] + "\n"
            desc_tmp = f.describe()
            res += str(desc_tmp[:10]) +"..." + str(desc_tmp[-10:]) +"\n "
        return res

