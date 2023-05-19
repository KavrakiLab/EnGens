'''
CrystalUtils - some cool utils for crystal structures workflow and analysis
'''
import os
import shutil
from typing import List
from enum import Enum
from os import path, system
from pathlib import Path
from xmlrpc.client import boolean

from Bio.PDB import PDBParser
from tqdm import tqdm
import mdtraj
import subprocess
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import SeqIO

class FastaSequence(object):

    def __init__(self, name: str, sequnce: List[str]) -> None:
        self.name = name
        self.sequence = sequnce

class CrystalUtils(object):

    # sequence utils
    d3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
            'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
            'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
            'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M', 'YCM':'C'}
    d1to3 = {v: k for k, v in d3to1.items()}
    
    def __init__(self, 
                pdb_codes: List[str] = None, 
                file_names: List[str] = None, 
                dst_folder: str = ".") -> None:

        # properties
        self.pdb_codes = pdb_codes
        self.file_names = file_names
        self.dst_folder = dst_folder
        self.dst_structure = Path(path.join(dst_folder, "structure_output"))
        self.dst_structure.mkdir(parents=True, exist_ok=True)
        self.dst_sequence = Path(path.join(dst_folder, "sequence_output"))
        self.dst_sequence.mkdir(parents=True, exist_ok=True)
        self.sequence_fasta = None
        self.msa_fasta = None

        # fetch files if only pdb codes were provided
        fetch_files = (file_names is None)
        if fetch_files: 
            self.file_names = []
            for pdb_code in tqdm(self.pdb_codes):
                if not len(pdb_code) == 4:
                    raise Exception("Found a pdb code that is not a 4 letter string: {}".format(pdb_code))
                if fetch_files:
                    file_res = self.fetchPDBFile(pdb_code, self.dst_structure)
                    self.file_names.append(file_res)
        else:
            # copy the files
            file_names_new = []
            for file_name in tqdm(self.file_names, desc="Copying files"):
                file_name_head, file_name_tail =  os.path.split(file_name)
                new_fileloc = path.join(self.dst_structure, file_name_tail)
                shutil.copy(file_name, new_fileloc)
                file_names_new.append(new_fileloc)
            self.file_names = file_names_new
            if pdb_codes is None or len(pdb_codes) == 0:
                self.pdb_codes = self.file_names 

    def extract_protein_sequence(self):
        # -----------------STEP 2 - extract protein sequences------------------------ #
        # Fasta file name - place to store the sequence
        experiment_name = "sequences"
        ff_name = path.join(self.dst_sequence, experiment_name+".fasta")

        pdb_seq_len=0
        # Create fasta file
        with open(ff_name, 'w') as ff_file:
            for i, pdb_f in enumerate(self.file_names):
                pdb_id = self.pdb_codes[i]
                pdb_seq = CrystalUtils.getSequence_1C(pdb_f)
                seq_rel = max(pdb_seq, key=len)
                ff_file.write(">{}\n".format(pdb_id))
                ff_file.write("{}\n".format(seq_rel))

        self.sequence_fasta = ff_name
        return ff_name
    
    def performMSA(self):
        
        if self.sequence_fasta is None:
            self.sequence_fasta = self.extract_protein_sequence()
        print("Running MSA of the sequences")
        # Location of ClustalO executable
        clust_path = "/clustalo-1.2.4-Ubuntu-x86_64" 
        in_file = self.sequence_fasta
        out_file = path.join(self.dst_sequence, "sequence_aligned.fasta")
        align_file = out_file
        if path.exists(align_file): 
            print("Found {}".format(align_file))
            print("Deleting old msa {}".format(align_file))
            os.remove(align_file)
        clustalomega_cline = ClustalOmegaCommandline(infile=in_file, outfile=out_file, verbose=True, auto=True, cmd=clust_path)
        print(clustalomega_cline)
        subprocess.run(str(clustalomega_cline).split())
        self.msa_fasta = align_file
        
        # extract the common regions
        # load the msa
        fasta_seqs = self.read_fasta(self.msa_fasta)
        cont_regs, pdb_elements_cont, msa_nonempty_pos = self.extract_common_regions()
        msa_vis_loc = self.visualizeMSA(self.msa_fasta, cont_regs)
        print("Visualization od MSA and common regions: {}".format(msa_vis_loc))

        print("Extracting substructures")
        #extract common substructures - only backbone and residues
        self.extract_common_substructures(pdb_elements_cont, msa_nonempty_pos, True)
        self.extract_common_substructures(pdb_elements_cont, msa_nonempty_pos, False)
    
    def extract_common_substructures(self, 
                                    pdb_elements_cont, 
                                    msa_nonempty_pos, 
                                    backbone_only: boolean = False):

        single_frames = []
        suffix = "_bbstrip.pdb" if backbone_only else "_resstrip.pdb"
        if backbone_only: 
            suffix = "_bbstrip.pdb"
            output_traj = "bb_traj.xtc"
            desc = "Extracting common regions from each file (backbone)"
            print(desc)
        else:
            suffix = "resstrip.pdb"
            output_traj = "resstrip_traj.xtc"
            desc = "Extracting common regions from each file (full residues)"
            print(desc)

        for pdb_elem in tqdm(pdb_elements_cont, desc=desc):
            output_bb_file = pdb_elem["file_loc"][:-4]+suffix
            print(output_bb_file)
            if not path.exists(output_bb_file):
                # read file, take subset of residues, save backbone 
                tmp_top = mdtraj.load(pdb_elem["file_loc"]).top
                allowed_residues = [pdb_elem["residue_idx_map"][x] for x in msa_nonempty_pos] 
                selstr = " or ".join(["resid == "+str(t) for  t in allowed_residues])
                if backbone_only:
                    selstr = " backbone and ({})".format(selstr)
                else:
                    selstr = " ({})".format(selstr)
                atom_sel = tmp_top.select(selstr)
                bb_loaded = mdtraj.load(pdb_elem["file_loc"], atom_indices= atom_sel)
                single_frames.append(bb_loaded)
                bb_loaded.save(output_bb_file)
            else:
                print("File exists: {}".format(output_bb_file))
        if backbone_only:
            print("Converting to trajectory {}".format(output_traj))
            cmd = ["mdconvert", "-f", path.join(self.dst_structure, "*"+suffix), "-o", path.join(self.dst_structure, output_traj)]
            subprocess.run(cmd)
    
    
    def extract_common_regions(self):
        aligned_file = self.msa_fasta
        print(aligned_file)
        fasta_sequences = SeqIO.parse(open(aligned_file),'fasta')
        # list of non-empty msa positions
        msa_nonEmpty_flags = None
        # dictionary id - non-empy indexes
        nonempty_idx_maps = {}
        #bookkeeping - see what amino-acids are common to everyone and can be selected in analysis
        for seq_i, f_seq in enumerate(fasta_sequences):
            seq = str(f_seq.seq)
            print(f_seq)
            head, tail = os.path.split(self.file_names[seq_i])
            key = tail
            # initialize the flags
            if msa_nonEmpty_flags is None:
                msa_nonEmpty_flags = [True for i in seq]
            # init the position map
            nonempty_idx_map = []
            # iterate through sequence
            # update msa nonempty map and nonempty index map
            cnt = 0
            for i, aa in enumerate(seq):
                if "-" in aa:
                    msa_nonEmpty_flags[i] = False
                    nonempty_idx_map.append(-1)
                else:
                    nonempty_idx_map.append(cnt)
                    cnt+=1
            nonempty_idx_maps[key] = nonempty_idx_map

        msa_nonEmpty_pos = [i for i, flag in enumerate(msa_nonEmpty_flags) if flag == True]
        # important - only these residues (global index) can besubprocess.run(cmd)
        # Show continuous regions
        cont_regs = []
        cont_reg = []
        for i, elem in enumerate(msa_nonEmpty_pos):
            if i == len(msa_nonEmpty_pos)-1:
                cont_reg.append(elem)
                cont_regs.append(cont_reg)
                break
            if msa_nonEmpty_pos[i+1]==msa_nonEmpty_pos[i]+1:
                cont_reg.append(elem)
            else: 
                if len(cont_reg)>0:
                    cont_reg.append(elem)
                    cont_regs.append(cont_reg)
                    cont_reg = []
                    print("Continuous region #{} found starting in AA range {}-{}".format(len(cont_regs), cont_regs[-1][0],
                                                                                        cont_regs[-1][-1]))

        # Save all this info in the structure with:
        # pdb_id
        # file location
        # list of allowed to select residues (0-indexed per file) (based on MSA)
        pdb_elements = []
        for i, pdb_id in enumerate(self.pdb_codes):
            head, tail = os.path.split(self.file_names[i])
            key = tail
            if not key in nonempty_idx_maps: continue
            pdb_struct = {}
            pdb_struct["pdb_id"] = pdb_id
            pdb_struct["file_loc"] = self.file_names[i]
            pdb_struct["residue_idx_map"] = nonempty_idx_maps[key]
            pdb_elements.append(pdb_struct)
        return cont_regs, pdb_elements, msa_nonEmpty_pos

    @staticmethod 
    def read_fasta(fasta_file:str) -> List[FastaSequence]:
        # Read in the resulting MSA
        mulseq_algn = []
        names = []
        curr_line = None
        fasta_seqs = []
        with open(fasta_file, "r") as file:
            for line in file.readlines():
                if line[0]==">":
                    if not curr_line is None:  
                        curr_line = curr_line.replace("\n", "")
                        seq = curr_line
                        mulseq_algn.append(seq)
                        fasta_seqs.append(FastaSequence(names[-1], seq))
                        name = line[1:-1]
                        names.append(name)
                    else:
                        name = line[1:-1]
                        names.append(name)
                    curr_line = ""
                else:
                    curr_line+=line
            
            fasta_seqs.append(FastaSequence(names[-1], seq))
        return fasta_seqs

    @staticmethod
    def getSequence_1C(pdb_file: str) -> List[str]:
        # Just an example input pdb
        record = pdb_file
        # run parser
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('struct', record)    
        # iterate each model, chain, and residue
        # printing out the sequence for each chain
        seqs = []
        for model in structure:
            for chain in model:
                seq = ""
                for residue in chain:
                    if not residue.resname in CrystalUtils.d3to1:
                        break
                    seq+=CrystalUtils.d3to1[residue.resname]
                seqs.append(seq)
        return seqs

    @staticmethod
    def fetchPDBFile(pdb_code:str, dst_folder: str = "."):
        pdb_file = path.join(dst_folder, pdb_code+".pdb")
        pdb_file_fixed = path.join(dst_folder, pdb_code+"_fixed.pdb")
        print("Fetching and fixing: {}".format(pdb_code))
        if path.exists(pdb_file_fixed):
            print("Already found")
            return pdb_file_fixed
        # fetch clean files
        pdbtools_command = "pdb_fetch {} | pdb_delhetatm | pdb_tidy > {}".format(pdb_code,pdb_file)
        system(pdbtools_command)
        # fix them still with PDBFixer
        pdbfixer_command = "/miniconda3/bin/python /usr/local/bin/pdbfix.py -i {} -o {}".format(pdb_file,pdb_file_fixed)
        system(pdbfixer_command)
        return pdb_file_fixed

    @staticmethod
    def visualizeMSA(fasta_file, regions=[]):

        # Read in the resulting MSA
        mulseq_algn = []
        names = []
        curr_line = None
        fasta_seqs = []
        with open(fasta_file, "r") as file:
            for line in file.readlines():
                if line[0]==">":
                    if not curr_line is None:  
                        curr_line = curr_line.replace("\n", "")
                        seq = curr_line
                        mulseq_algn.append(seq)
                        fasta_seqs.append(FastaSequence(names[-1], seq))
                        name = line[1:-1]
                        names.append(name)
                    else:
                        name = line[1:-1]
                        names.append(name)
                    curr_line = ""
                else:
                    curr_line+=line
            mulseq_algn.append(seq)
            
        # generate JS txt filler
        #sequences
        js_txt = "const seqs = [ "
        for i, pdb_id in enumerate(names):
            js_txt += "{ name : \""+ pdb_id +"\", sequence: \""+mulseq_algn[i]+"\"}, \n"
        js_txt = js_txt[:-3]+"]; \n"
        js_txt += "const sequence = \""+mulseq_algn[0]+"\"; \n"
        
        # extracted regions
        output_name = path.splitext(fasta_file)[0]+".html"
        if len(regions) > 0:
            output_name = path.splitext(fasta_file)[0]+"_regions.html"
            js_txt+="const features=[ {\"accession\": \"extractedRegions\",  \"color\": \"#342ea2\", \"locations\": [  {\"fragments\": [ \n"
            for i, region in enumerate(regions):
                js_txt += "{ \
                \"start\": "+str(region[0]+1)+", \
                \"end\": "+str(region[-1]+1)+"}, \n"
            js_txt = js_txt[:-3]+" ] } ] } ];\n"
        else:
            js_txt+="const features=[];\n"

        seq_length = len(mulseq_algn[0])
        with open("/usr/local/template_msa.html", "rt") as fin:
            with open(output_name, "wt") as fout:
                for line in fin:
                    new_line = line.replace('$VARIABLES_PLACEHOLDER$', js_txt)
                    new_line = new_line.replace('$SL$', "\"{}\"".format(seq_length))
                    fout.write(new_line)
        return output_name

# -------------------------------------------------- #


