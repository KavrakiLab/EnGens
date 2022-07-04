
from math import ceil
import pytraj as pt
import nglview as ngl
import mdtraj
import pyemma.coordinates 
import pyemma.coordinates.data
import os
import shutil
import numpy as np
import tqdm
import importlib
import nglview
from Bio import PDB
from subprocess import call

'''
EnGen - the main class for ensemble generation and analysis 
'''

class EnGen(object):

    def __init__(self, trajectory, topology, 
                topology_select=None,
                cryst_pdb_list = False,
                file_names = None, 
                align=False, 
                chunk_size=1000):
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
        cryst_pdb_list: boolean, default = False
            is your input a list of crystal structure files??
        file_names: list, default = None
            if your input is the list of PDB files, provide the list here
        align: boolean, default = False
            specify if you wish to align the given trajectory before analyzing
        chunk_size: int, default = 1000
            specify the chunk size of trajectory that is loaded into memory at once
            (decrease this number when having memory issues)
        """

        # perform the alignment if specified

        self.traj = trajectory
        self.traj_name = trajectory
        self.full_traj_name = trajectory
        self.ref = topology
        self.full_ref = topology
        self.refrestraint = topology_select
        self._selection_indices = None
        self.chunk_size = chunk_size
        self.mdtrajref = mdtraj.load(self.ref).topology
        self.crystal_flag = cryst_pdb_list
        self.pdb_files = file_names
        self.pdb_list = None
        self.mdtraj_list = []

        if align or cryst_pdb_list:
            traj_new = trajectory[:trajectory.rfind(".")]+"-aligned.xtc"
            self.align_trajectory(traj_new)
            self.traj = traj_new
            self.traj_name = traj_new
            self.full_traj_name = traj_new

        if align and cryst_pdb_list:
            self.align_pdb_files()


        if self.crystal_flag and self.refrestraint is None:
            self.pdb_list = []
            self.mdtraj_list = []
            for elem in tqdm.tqdm(self.pdb_files, "Loading files (might take a while)"):
                elem_frame = mdtraj.load(elem)
                elem_frame_dst = elem[:-4]+"_tmp.xtc"
                elem_frame.save(elem_frame_dst)
                self.mdtraj_list.append(elem_frame)
                pyemma_frame = pyemma.coordinates.source(elem_frame_dst, top=elem)
                self.pdb_list.append(pyemma_frame)


        # load a restrained trajectory based on a selection/list/None
        if isinstance(self.refrestraint, str):
            tmp_top = mdtraj.load(self.ref).topology
            self._selection_indices = tmp_top.select(topology_select)
            self.select_atoms_trajectory(self._selection_indices)

        elif isinstance(self.refrestraint, list):
            self._selection_indices = self.refrestraint
            self.select_atoms_trajectory(self._selection_indices)


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
    

    def align_pdb_files(self):

        parser = PDB.PDBParser(QUIET = True)
        # first structure as reference structure
        ref_structure = parser.get_structure("tmp_ref", self.pdb_files[0])
        ref_atoms = [] # only align C- alpha
        # Iterate of all chains in the model in order to find all residues
        for ref_chain in ref_structure[0]:
            # Iterate of all residues in each model in order to find proper atoms
            for ref_res in ref_chain:
                ref_atoms.append(ref_res['CA'])

        for elem in tqdm.tqdm(self.pdb_files, "Aligning pdb files (might take a while)"):
            sample_structure = parser.get_structure("tmp_sample", elem)
            sample_model = sample_structure[0]
            sample_atoms = [] # only Calpha
            for sample_chain in sample_model:
                for sample_res in sample_chain:
                    sample_atoms.append(sample_res['CA'])
            # Now we initiate the superimposer:
            super_imposer = PDB.Superimposer()
            super_imposer.set_atoms(ref_atoms, sample_atoms)
            super_imposer.apply(sample_model.get_atoms())
            # Save the aligned version of 1UBQ.pdb
            io = PDB.PDBIO()
            io.set_structure(sample_structure) 
            io.save(elem[:-4]+"_algn.pdb")
        
        self.pdb_files = [ f[:-4]+"_algn.pdb" for f in self.pdb_files]

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

    def select_atoms_pdb_list(self, selected_atoms):

        self.pdb_list = []
        old_files = self.pdb_files
        self.pdb_files = []
        self.mdtraj_list = []
        for elem in tqdm.tqdm(old_files, "Loading PDB files with atom selection"):
            elem_frame = mdtraj.load(elem, atom_indices= selected_atoms)
            elem_frame_dst = elem[:-4]+"_tmp.xtc"
            elem_frame.save(elem_frame_dst)
            pyemma_frame = pyemma.coordinates.source(elem_frame_dst, top=elem_frame)
            self.mdtraj_list.append(elem_frame)
            self.pdb_list.append(pyemma_frame)
            elem_frame_dst = elem[:-4]+"_tmp.pdb"
            elem_frame.save(elem_frame_dst)
            self.pdb_files.append(elem_frame_dst)
    
    def select_atoms_trajectory(self, selected_atoms):
        
        if self.crystal_flag:
            self.select_atoms_pdb_list(selected_atoms)
        tmp_traj = mdtraj.iterload(self.traj, top=self.ref, atom_indices=selected_atoms)
         # see trajectory length
        tmp_pe = pyemma.coordinates.source(self.traj, top=self.ref)
        traj_len = tmp_pe.trajectory_length(0)
        file_names = []
        n_iter = ceil(traj_len/self.chunk_size)
        saved_ref=False
        tmp_dir ="./tmp_files"
        if not os.path.exists(tmp_dir):
            os.makedirs(tmp_dir)
        else:
            shutil.rmtree(tmp_dir)
            os.makedirs(tmp_dir)

        i=0
        # save chunks
        for chunk in tqdm.tqdm(tmp_traj, desc = "Making the selection... "):

            if not saved_ref:
                tmp_pdb_name = self.traj[:self.traj.rfind(".")]+"-engen-selected.pdb"
                chunk[0].save(tmp_pdb_name)
                self.ref = tmp_pdb_name
                self.mdtrajref = mdtraj.load(tmp_pdb_name).topology
                self.ref = tmp_pdb_name
                saved_ref=True
            # save a temporary chunk
            file_name = os.path.join(tmp_dir, "engen-tmp_chunk"+str(i)+".xtc")
            chunk.save(file_name)
            file_names.append(file_name)
            i+=1 
        
        tmp_name = self.traj[:self.traj.rfind(".")]+"-engen-selected.xtc"
        
        if os.path.exists(tmp_name):
            os.remove(tmp_name)

        ret_val = call(["mdconvert -o {} {}".format(tmp_name, " ".join(file_names))] , shell=True)
        if not ret_val == 0:
            raise(Exception("Error making selection."))

        self.traj = tmp_name
        self.traj_name = tmp_name

        for i in tqdm.tqdm(range(len(file_names)), desc="Cleaning files..."):
            elem = file_names[i]
            os.remove(elem)

        if len(os.listdir(tmp_dir))==0:
            os.rmdir(tmp_dir)

    def align_trajectory(self, output_name):

        traj_name = self.traj_name 
        ref_name = self.ref
        # create the iterator
        mdl = mdtraj.iterload(traj_name, top=ref_name, chunk=self.chunk_size)
        # see trajectory length
        tmp_pe = pyemma.coordinates.source(traj_name, top=ref_name)
        traj_len = tmp_pe.trajectory_length(0)
        first_frame = None
        file_names = []
        n_iter = ceil(traj_len/self.chunk_size)
        tmp_dir ="./tmp_files"
        if not os.path.exists(tmp_dir):
            os.makedirs(tmp_dir)
        else:
            shutil.rmtree(tmp_dir)
            os.makedirs(tmp_dir)
        # align and save chunks
        for i in tqdm.tqdm(range(n_iter), desc = "Aligning trajectory: "):
            chunk = next(mdl)
            if first_frame is None:
                first_frame = chunk[0]
            chunk_s = chunk.superpose(first_frame)
            # save a temporary aligned chunk
            file_name = os.path.join(tmp_dir, "engen-tmp_chunk"+str(i)+".xtc")
            chunk_s.save(file_name)
            file_names.append(file_name)

        if os.path.exists(output_name):
            os.remove(output_name)

        ret_val = call(["mdconvert -o {} {}".format(output_name, " ".join(file_names))] , shell=True)
        if not ret_val == 0:
            raise(Exception("Error making alignment."))

        for i in tqdm.tqdm(range(len(file_names)), desc="Cleaning files..."):
            elem = file_names[i]
            os.remove(elem)

        if len(os.listdir(tmp_dir))==0:
            os.rmdir(tmp_dir)


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
        if self.crystal_flag:
            pyemma_feat = []
            for elem in tqdm.tqdm(self.pdb_files, "Adding featurizers per PDB file.."):
                name = ""
                tmp_top = pyemma.coordinates.featurizer(elem)
                for key, params in feats.items():
                    name+=key[len("add_"):]
                    name+="&"
                    func =  getattr(tmp_top, key)
                    func(**params)
                name = name[:-1]
                pyemma_feat.append(tmp_top)
            self.featurizer_names.append(name)
            self.featurizers.append(pyemma_feat)
        else:
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

        if not self.crystal_flag:
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

        else:
            
            n_residues = min(mdtraj.load(file).top.n_residues for file in self.pdb_files)
            all_resi = [i for i in range(n_residues)]
            default_feat1 = {
                "add_residue_COM": {"residue_indices":all_resi ,"scheme":'all'}
            }
            default_feats = [default_feat1]

        for feat_dict in default_feats:
            self.add_featurizer(feat_dict)
        
    def apply_featurizations(self):
        self.data = []
        for f in self.featurizers:
            if not self.crystal_flag:
                self.data+=pyemma.coordinates.source(self.traj, features=f)
            else:
                per_pdb_data = None
                for i, elem in tqdm.tqdm(enumerate(self.mdtraj_list), "Applying featurization per PDB file"):
                    tmp_name = self.pdb_files[i][:-4]+".xtc"
                    elem.save(tmp_name)
                    data = pyemma.coordinates.load(tmp_name, features=f[i])
                    if per_pdb_data is None: 
                        per_pdb_data = data
                    else:
                       per_pdb_data = np.vstack((per_pdb_data, data)) 
                self.data.append((0, per_pdb_data))


    def clean(self):
        if self.crystal_flag == True:
            # free this space
            self.mdtraj_list = None
            self.pdb_list = None

    def describe_featurizers(self):
        res = ""
        if not self.crystal_flag:
            for i, f in enumerate(self.featurizers):
                res += "Featurizer no. "+str(i)+":\n "
                res += self.featurizer_names[i] + "\n"
                desc_tmp = f.describe()
                res += str(desc_tmp[:10]) +"..." + str(desc_tmp[-10:]) +"\n "
        else:
            for i, f in enumerate(self.featurizers):
                res += "Featurizer no. (residues may differ from file to file) "+str(i)+":\n "
                res += self.featurizer_names[i] + "\n"
                desc_tmp = f[i].describe()
                res += str(desc_tmp[:10]) +"..." + str(desc_tmp[-10:]) +"\n "

        return res

js_script = """
var x = document.nglview.stage.getRepresentationsByName("selection");
var stickRepr = x['list'][0];
var rules = JSON.stringify(stickRepr.repr.selection.selection.rules);
console.log("Hello");
console.log(rules);
var command = "selection = '" + rules + "'";
IPython.notebook.kernel.execute(command);
IPython.notebook.kernel.execute("selection = json.loads(selection)");
"""

def get_selstring(selection):
    chains = []
    residues = []
    for elem in selection:
        for rule in elem['rules']:
            for key, value in rule.items():
                if key == "chainname": chains.append(value)
                if key == "resno": residues.append(value)

    sel_string = ""
    for residue in residues:
        sel_string+= "residue=="+str(residue) + " or "
    sel_string = sel_string[:-len(" or ")]
    return sel_string

def select_residues_nglview(top_loc):
    nglwidget = nglview.show_structure_file(top_loc)
    nglwidget.clear_representations()
    nglwidget.add_cartoon(colorScheme="residueindex")
    nglwidget.add_ball_and_stick(color="red", selection="0", name="selection")
    nglwidget.gui_style = 'ngl'
    nglwidget._execute_js_code("document.nglview = this;")
    nglwidget._execute_js_code(
    """
        var stickSel = ""
        var x = this.stage.getRepresentationsByName("selection")
        var stickRepr = x['list'][0]

        var f1 = function (pickingProxy) {
        if (stickRepr.repr.selection.selection.rules[0] && stickRepr.repr.selection.selection.rules[0].keyword == 20) {
         stickSel = ""
        }
        if (pickingProxy && pickingProxy.ctrlKey && (pickingProxy.atom || pickingProxy.bond)){
            console.log("CTRL")
            console.log(pickingProxy)
            var atom = pickingProxy.atom || pickingProxy.closestBondAtom;
            var residue = atom.residue
            var curSel = String(residue.resno)+':'+residue.chainname+' or '

            console.log(curSel)

            var isSel = stickSel.search(curSel)
            if (isSel == -1) {
                // Append to selection
                stickSel += curSel
            }
            console.log(stickSel);
            stickRepr.setSelection(stickSel)

        }

        if (pickingProxy && pickingProxy.shiftKey && (pickingProxy.atom || pickingProxy.bond)){
            console.log("SHIFT")
            console.log(pickingProxy)
            var atom = pickingProxy.atom || pickingProxy.closestBondAtom;
            var residue = atom.residue
            var curSel = String(residue.resno)+':'+residue.chainname+' or '
            console.log(curSel)
            console.log(stickSel)
            var isSel = stickSel.search(curSel)
            if (isSel != -1)  {
                // Remove from selection
                stickSel = stickSel.replace(curSel, "")
            }
            console.log(stickSel);
            if(stickSel.length == 0) {
                stickRepr.setSelection("none")
            }
            else{
            stickRepr.setSelection(stickSel)
            }
        }

        }
        this.stage.signals.hovered.add(f1)
    """
    )
    return nglwidget