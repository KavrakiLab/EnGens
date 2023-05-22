from pathlib import Path
import Bio
from Bio.PDB import PDBList
import shutil
from Bio.PDB import PDBParser
import mdtraj
from Bio.Align.Applications import ClustalOmegaCommandline
import subprocess
from Bio import SeqIO
import matplotlib.pyplot as plt
from os import path, system
from engens.core.CrystalUtils import *
import requests
import pandas as pd
import numpy as np
import os
import sys
from tqdm import tqdm
from pathlib import Path
from os import path
from engens.core.cryst_utils import *
from enum import Enum
import math

INPUT_TYPE = Enum("INPUT_TYPE", ["UNIPROT", \
                                "PDB_CODES", \
                                "PDB_FILES"])

class PrepStatic:
    
    def __init__(self, 
            pdb_codes: List[str] = None, 
            file_names: List[str] = None,
            uniprot_ids: List[str] = None, #- TODO
            dst_folder: str = ".") -> None:

        # properties
        self.input_type = None
        self.uniprot_ids = uniprot_ids
        self.pdb_codes = pdb_codes
        self.file_names = file_names
        self.dst_folder = dst_folder
        self.dst_structure = Path(path.join(dst_folder, "structure_output"))
        self.dst_structure.mkdir(parents=True, exist_ok=True)
        self.dst_sequence = Path(path.join(dst_folder, "sequence_output"))
        self.dst_sequence.mkdir(parents=True, exist_ok=True)
        self.uniprot_metadata = None
        self.pdb_metadata = None
        self.domains = []
        self.mTM_align_per_domain = {}
        self.domain_metadata = None
        self.final_traj_file = None
        self.final_pdb_files = None

        #---------------STEP 1 - LOAD THE DATA------------------#

        # if input is PDB codes
        if self.uniprot_ids is not None:
            self.input_type = INPUT_TYPE.UNIPROT
            print("======================================================")
            print("STEP 0 - Querying UniProt DB for the IDs")
            print("======================================================")
            uniprot_details, pdb_details = uniprot_get_details(self.uniprot_ids, True)
            print()
            print("Found structures:")
            
            print(tabulate(pdb_details,
                          headers = pdb_details.columns))
            print("Attention - this step requires user input! (press enter to continue)")
            input()   
            selected_pdbs = pick_PDBs(pdb_details)
            self.pdb_codes = list(selected_pdbs.pdb_id.unique())
            
        # if input is PDB codes
        if self.pdb_codes is not None:
        #---------------STEP 1A - LOAD THE DATA FROM PDB CODES------------------#
            self.input_type = INPUT_TYPE.PDB_CODES
            # 1A.1 download them from pdb renum and fix them
            file_names = {}
            problem_download_pdbs = []
            print("======================================================")
            print("STEP 1 - Downloading renumbered pdbs and fixing files")
            print("======================================================")
            for pdb_id in tqdm(self.pdb_codes):
                res_file = get_pdb_bioassembly_pdbrenum(pdb_id.lower(), self.dst_structure)
                if res_file is None:
                    print("Could not download {} - discarding pdb".format(pdb_id))
                    problem_download_pdbs.append(pdb_id)
                else:
                    file_name = res_file[:-5]+"_fixed.pdb"
                    fix_pdb_files(res_file, file_name)
                    file_names[pdb_id] = file_name
            self.pdb_codes = [x for x in self.pdb_codes if not x in problem_download_pdbs]
            print("======================================================")
            print("Successfull download and fixing of files. Location: {}".format(self.dst_structure))
            print("======================================================")
            # 1A.2 get their UniProt and PDB metadata
            print("======================================================")
            print("STEP 2 - Downloading matedata associated with given codes")
            print("======================================================")
            print()
            print("Fetching the pdb and uniprot metadata for codes")
            print()
            
            self.pdb_metadata, self.uniprot_metadata = get_pdb_metadata(self.pdb_codes)
            print()
            print("PDB code associated metadata:")
            print(tabulate(self.pdb_metadata,
                          headers = self.pdb_metadata.columns))
            self.pdb_metadata["file_loc"] = self.pdb_metadata["pdb_id"].apply(lambda x: file_names[x])
            print()
            print("UNIPROT metadata:")
            print(tabulate(self.uniprot_metadata[["accession_id", "id", "full_name"]], \
                           headers = ["accession_id", "id", "full_name"]))
            print()
            print("======================================================")
            print("STEP 3 - Define domains and map UNIPROT accessions to the domains")
            print("======================================================")
            print("Attention - this step requires user input! (press enter to continue)")
            input()
            # 1A.3 select main domains/chains and their mapping to uniprot ids within each pdb 
            self.uniprot_metadata, self.pdb_metadata, domain_metadata = define_domains(self.uniprot_metadata, self.pdb_metadata)
            self.domains = list(domain_metadata.domain_name.unique())
            self.domain_metadata = domain_metadata
            
            print()
            print("UNIPROT metadata mapped to domains:")
            print(tabulate(self.uniprot_metadata[["domain", "accession_id", "id", "full_name"]], \
                           headers = ["domain", "acc_id", "id", "full_name"]))
            
            print()
            print("DOMAIN metadata:")
            print(tabulate(self.domain_metadata, \
                           headers = self.domain_metadata.columns))
            
            print()
            print("======================================================")
            print("STEP 4 - Defining final selection for processing")
            print("======================================================")
            print()
            print("Attention - this step requires user input! (press enter to continue)")
            input()
            n_domains = self.domain_metadata.domain_name.unique().shape[0]
            print("Your current inputs contain {} domains with the following associated uniprots".format(n_domains))

            tmp_df = pd.merge(self.domain_metadata,self.uniprot_metadata, \
                              left_on = "uniprot_id", right_on="accession_id")\
                            [["domain_name", "uniprot_id", "full_name"]]
            columns = ["domain_name", "uniprot_id", "accession_name"]
            print(tabulate(tmp_df, columns))

            print()
            print("Please choose the domain-accession pairs you want to consider for your main analysis ")

            while True:
                in_sele = input("Input indices of domain-uniprot metadata separated by space (e.g., '0 3 5' )")
                try: 
                    in_list = [int(x) for x in in_sele.split(" ")]
                    selected_unis = tmp_df.iloc[in_list]
                except:
                    print("Invalid input, try again")
                    continue
                print("Selected: ")
                print(tabulate(selected_unis[["domain_name", "uniprot_id", "full_name"]], \
                               headers = columns))

                selected_dom_acc_ids = list(selected_unis.apply(lambda x: x.domain_name+"-"+x.uniprot_id, axis=1))
                print(selected_dom_acc_ids)
                break

            print()
            print("PDB codes satisfying the above selection (containing all domains with any of the related uniprot accessions)")
            self.selected_unis = selected_unis
            selected_pdbs = find_pdb_ids_dom_uni_qualifying(self.pdb_metadata, selected_unis)
            print()
            print(selected_pdbs.pdb_id.unique())
            print("With following metadata:")
            print(tabulate(selected_pdbs))
            print()
            print("Discarding the following entries (that do not contain selected domains:")
            discarded_pdbs = set(list(self.pdb_metadata.pdb_id.unique()))-set(list(selected_pdbs.pdb_id.unique()))
            print(discarded_pdbs)
            self.pdb_metadata = self.pdb_metadata[self.pdb_metadata["pdb_id"].isin(list(selected_pdbs.pdb_id.unique()))]
            #self.pdb_metadata = self.pdb_metadata.merge(selected_pdbs, left_on = ["pdb_id", "domain", "accession"],\
            #                                                     right_on = ["pdb_id", "domain", "uniprot_id"], how="inner")
            self.domain_metadata = self.domain_metadata.merge(selected_unis, on = ["domain_name", "uniprot_id"],\
                                                                        how="inner")
            self.pdb_metadata = self.pdb_metadata[~self.pdb_metadata.domain.isna()]
            self.pdb_metadata = self.pdb_metadata.drop_duplicates(subset=["pdb_id", "accession", "domain"])
            print()
            print("======================================================")
            print("STEP 5- Extracting coordinates associated with given domains (per-accession)")
            print("======================================================")
            # 1A.4 Extract domain atoms into separate clean pdb files
            self.domain_metadata["pdb_files_dir"] = self.domain_metadata.apply(lambda x: Path(path.join(self.dst_structure, \
                                                                          x.domain_name+"/"+x.uniprot_id)),\
                                                                  axis=1)
            self.domain_metadata["pdb_files_dir"].apply(lambda x: x.mkdir(parents=True, exist_ok=True))
            
            self.pdb_metadata["domain_file_loc"] = self.pdb_metadata.apply(lambda x: extract_chains_create_file(x.file_loc, \
                                                                                                                x.domain, \
                                                                                                                x.accession, \
                                                                                                                x.first_asym_id), \
                                                                                                                axis=1)
            
            self.pdb_metadata= self.pdb_metadata[~self.pdb_metadata["domain_file_loc"].isna()]
            self.pdb_metadata["domain_file_loc"].apply(lambda x: pdb_reorder_residues(x))
            
            print("Extracted per-accession domain coordinates")
            print("DOMAIN metadata:")
            print(tabulate(self.domain_metadata, \
                           headers = self.domain_metadata.columns))
            print("======================================================")
            print("STEP 6- Creating domain alignments (accoarding to uniprot accession numbering)")
            print("======================================================")
            # 1A.5 Align per domains(per uniprot id)
            
            self.domain_metadata["align_dir"] = self.domain_metadata.apply(lambda x: Path(path.join(self.dst_sequence, \
                                                                         x.domain_name+"/"+x.uniprot_id)),\
                                                              axis=1)
            self.domain_metadata["align_dir"].apply(lambda x: x.mkdir(parents=True, exist_ok=True))
            
            starts = []
            ends = []
            res_files = []
            for i, row in self.domain_metadata.iterrows():
                fasta_res_file = path.join(row.align_dir, "msa.fasta")
                start, end, _ = extract_msa_for_domain_signle_uniprot(row.domain_name, \
                                                                      row.uniprot_id, \
                                                                      fasta_res_file, \
                                                                      self.pdb_metadata)
                starts.append(start)
                ends.append(end)
                res_files.append(fasta_res_file)
                visualizeMSA(fasta_res_file)
                print("Check out the MSA visualized at {}.html".format(fasta_res_file[:-6]))
                
            self.domain_metadata["start"] = starts
            self.domain_metadata["end"] = ends
            self.domain_metadata["fasta_file"] = res_files
            
            print("Extracted per-accession domain alignments")
            #print("DOMAIN metadata:")
            #print(tabulate(self.domain_metadata, \
            #               headers = self.domain_metadata.columns))
            print("======================================================")
            print("STEP 7- Extracting per-accession domain MCS")
            print("======================================================")
            
            # 1A.6 Find MCS and MCS regions per domain-accession pairs
            self.domain_metadata["mcs"] = self.domain_metadata["fasta_file"].apply(lambda x: \
                                                                                   get_mcs_residues_from_fasta(x))
            self.domain_metadata["mcs_regions"] = self.domain_metadata["mcs"].apply(lambda x: \
                                                                                    [i for i in intervals_extract(x)])
            
            self.domain_metadata.apply(lambda x: visualizeMSA(x.fasta_file, x.mcs_regions),axis=1)
            
            # 1A.7 Extract MCS portions of the files
            self.domain_metadata.apply(lambda x: extract_mcs(x.pdb_files_dir, x.mcs, x.start), axis=1)
            # check if the extracted thing has a consistent backbone
            self.domain_metadata.apply(lambda x: check_atom_numbers_mcs_match(x.pdb_files_dir,), axis=1)
            
            
            print("Successfully extracted per-accession domain MCS!")
            print("======================================================")
            print("STEP 8 - Processing multi-accession domains (requires mTM-align)")
            print("======================================================")
            
            self.mcs_metadata = perform_mTM_align(self.domain_metadata, self.pdb_metadata, self.dst_structure)
            
            final_mcs_dst = path.join(self.dst_structure, "final_mcs")
            Path(final_mcs_dst).mkdir(parents=True, exist_ok=True)
            print("Successfully processed multi-accession domain MCS! (if present)")
            print("======================================================")
            print("STEP 9 - Combining the final MCS data")
            print("======================================================")
            # if no mTM-align is necessary, just copy the uniprot mcs files to final location
            iterrows_object = None
            if self.mcs_metadata is None:
                print()
                print("Not the case of multi-accession - MCS extracted based on accession residue numbers")
                iterrows_object = self.domain_metadata.iterrows()
            else:
                print()
                print("Multi-accession - MCS extracted based on accession mTM-align")
                iterrows_object = self.mcs_metadata.iterrows()
                
            res_files = {}
            res_pdb_codes = {}
            print()
            print("Step 9.1 - copy all MCS to {} for processing and extract backbones".format(final_mcs_dst))
            for i, row in iterrows_object:
                
                if self.mcs_metadata is None:
                    dom = row.domain_name
                else:
                    dom = row.domain
                
                if self.mcs_metadata is not None:
                    resn = row.mcs_resn
                    
                uni = row.uniprot_id
                
                res_files[dom+"-"+uni] = []
                pdbs_mtd = self.pdb_metadata[(self.pdb_metadata.domain == dom) & \
                                                 (self.pdb_metadata.accession == uni)]
                mcs_files = pdbs_mtd.domain_file_loc.apply(lambda x: x[:-4]+"_mcs.pdb")
                pdb_codes = list(pdbs_mtd.pdb_id)
                for pdb_code in pdb_codes:
                    if not pdb_code in res_pdb_codes:
                        res_pdb_codes[pdb_code] = [dom+"-"+uni]
                    else:
                        res_pdb_codes[pdb_code].append(dom+"-"+uni)

                # extract the backbones of the files 
                for file in mcs_files:
                    new_file = path.join(final_mcs_dst, dom+"_"+uni+"_"+Path(file).stem+".pdb")
                    os.system("cp {} {}".format(file, new_file))
                    new_file_bb = path.join(final_mcs_dst, dom+"_"+uni+"_"+Path(file).stem+"_bb.pdb")
                    #extract backbone

                    tmp = pd.read_csv(path.join(file), sep="\n", header=None)
                    tmp["res_n"] = tmp[0].apply(get_residue_number)
                    tmp["res_id_3"] = tmp[0].apply(get_residue_id)
                    tmp["atom"] = tmp[0].apply(lambda x: x.split()[2])
                    tmp = tmp[tmp["atom"].isin(["N", "CA", "C", "O"])]
                    if self.mcs_metadata is not None:
                        tmp = tmp[tmp["res_n"].isin(resn)]
                    tmp[0].to_csv(new_file_bb, sep="\n", index=False, header=False)
                    res_files[dom+"-"+uni].append(new_file_bb)

            print("Successfully copied all cleaned MCS domains to {}!".format(final_mcs_dst))
            print()
            print("Step 9.2 - combine multiple domains from same pdb id to single files")
            # combine domains into files 
            final_files = []
            for pdb_id_key in res_pdb_codes:
                print(pdb_id_key)
                dom_uni_entries = res_pdb_codes[pdb_id_key]
                print(dom_uni_entries)
                files_to_concat = []
                for dom_uni in dom_uni_entries:
                    tmp_dom_uni_files = res_files[dom_uni]
                    for file in tmp_dom_uni_files:
                        if pdb_id_key in file or pdb_id_key.lower() in file:
                            files_to_concat.append(file)
                final_file = path.join(final_mcs_dst, pdb_id_key+"_mcs_bb.pdb")
                print(files_to_concat)
                os.system("cat {} > {}".format(" ".join(files_to_concat), final_file))
                final_files.append(final_file)
            print("Successfully combined all cleaned MCS domains to complex PDBs in {}!".format(final_mcs_dst))
            # combine into the trajectory

            print()
            print("Step 9.3 - combine multiple files to single EnGens input")
            self.final_traj_file = path.join(final_mcs_dst, "final_bb_traj.xtc")
            self.final_pdb_files = final_files
            os.system("rm {}".format(self.final_traj_file))
            os.system("mdconvert {} -o {}".format(" ".join(final_files), self.final_traj_file))
            print("Successfully combined all cleaned MCS domains to complex PDBs in {}!".format(final_mcs_dst))

            print("======================================================")
            print("Successfully completed PDB pre-processing")
            print("======================================================")
                    
            # extract the MCS and all related data
            # TODO
            # 1 - extract all components and create a mock trajectory
            
            # 2 - process the MSA to showcase the final MCS
            
        
        # if input are PDB files just fix them
        # they will be treated in bulk as a single domain
        elif self.file_names is not None:
            # TODO
            """
            self.input_type = INPUT_TYPE.PDB_FILES
            print("Inspecting and fixing files")
            for file_name in file_names:
                file_name_stem = Path(file_name).stem
                file_name_dir = Path(file_name).parent
                fix_pdb_files(file_name, path.join(file_name_dir, file_name_stem+"_fixed.pdb"))
            print("No metadata available for files")
            
            #mTM-align is necessary for MCS as no metadata is available
            self.domains = ["DOMAIN0"]
            self.mTM_align_per_domain["DOMAIN0"] = True
            #self.domain_map = domain_map
            """
        # if input are UniProt Ids, fetch the corresponding PDBs
        ##TODO         