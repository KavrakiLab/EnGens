import pandas as pd
import re
from engens.core.CrystalUtils import *
import Bio.PDB.PDBParser
import numpy as np

#converts the string with gaps into list of positions of aa-s
def alignment_string_to_numbers(alignment):
    num_array = []
    cnt = 0
    for aa in alignment:
        if aa.isalpha():
            num_array.append(cnt)
            cnt+=1
        else:
            cnt+=1
            num_array.append(-1)
    return num_array

import itertools
  
def intervals_extract(iterable):
    iterable = sorted(set(iterable))
    for key, group in itertools.groupby(enumerate(iterable),
    lambda t: t[1] - t[0]):
        group = list(group)
        yield [group[0][1], group[-1][1]]

        
from Bio import Align

def multi_uniprot_alignment_pairwise(meta_uniprot_mapping, uniprot_details):

    meta_uniprot_data = {
        "query_id":[],
        "target_id":[],
        "query_seq_aligned":[],
        "target_seq_aligned":[],
        "mapping":[],
        "reverse_mapping":[]
    }

    global_sequence_mappings = {}
    global_sequence_mappings_reverse = {}

    for query_id, target_id in meta_uniprot_mapping.items():
        query_sequence = uniprot_details[uniprot_details.accession_id==query_id]["seq"].iloc[0]
        target_sequence = uniprot_details[uniprot_details.accession_id==target_id]["seq"].iloc[0]
        meta_uniprot_data["query_id"].append(query_id)
        meta_uniprot_data["target_id"].append(target_id)
        
        if query_id == target_id:
            query_ids = list(range(len(query_sequence)))
            global_sequence_mappings[query_id] = dict(zip(query_ids, query_ids))
            global_sequence_mappings_reverse[query_id] = dict(zip(query_ids, query_ids))
            meta_uniprot_data["query_seq_aligned"].append(query_sequence)
            meta_uniprot_data["target_seq_aligned"].append(query_sequence)
            meta_uniprot_data["mapping"].append(dict(zip(query_ids, query_ids)))
            meta_uniprot_data["reverse_mapping"].append(dict(zip(query_ids, query_ids)))
        else:    
            aligner = Align.PairwiseAligner()
            aligner.mode = 'global'
            alignments = aligner.align(target_sequence, query_sequence)
            target_res = alignments[0].target
            query_res = alignments[0].query
            target_array = alignment_string_to_numbers(target_res)
            query_array = alignment_string_to_numbers(query_res)
            all_mcs_positions = list(range(len(query_aa_position_array)))
            global_sequence_mappings[query_id] = dict(zip(query_array, all_mcs_positions))
            global_sequence_mappings_reverse[query_id] = dict(zip(all_mcs_positions, query_array))
            meta_uniprot_data["query_seq_aligned"].append(query_res)
            meta_uniprot_data["target_seq_aligned"].append(target_res)
            meta_uniprot_data["mapping"].append(dict(zip(query_array, all_mcs_positions)))
            meta_uniprot_data["reverse_mapping"].append(dict(zip(all_mcs_positions, query_array)))
    return pd.DataFrame(meta_uniprot_data)

def multi_uniprot_alignment_from_files(file_name, meta_uniprots, meta_uniprot_mapping):
    
    seq_aligned = {}
    for meta_uniprot in meta_uniprots:
        fname = file_name.format(meta_uniprot)
        print(fname)
        name = None
        line_name = False
        with open(fname, 'r') as fasta_file:
            lines = fasta_file.readlines()
            for line in lines:
                if line[0] == ">":
                    if name is not None:
                        seq_aligned[name] = "".join(seq_aligned[name]).replace("\n", "")
                    name = line.split('|')[1]
                    line_name = True
                    print(name)
                else:
                    line_name = False 

                if not line_name:
                    if name not in seq_aligned:
                        seq_aligned[name] = [line]
                    else:
                        seq_aligned[name].append(line)
                    line_name = False
        if name is not None:
            seq_aligned[name] = "".join(seq_aligned[name]).replace("\n", "")
        
    # Meta alignment 
    meta_uniprot_data = {
        "query_id":[],
        "target_id":[],
        "query_seq_aligned":[],
        "target_seq_aligned":[],
        "mapping":[],
        "reverse_mapping":[]
    }

    global_sequence_mappings = {}
    global_sequence_mappings_reverse = {}

    for query_id, target_id in meta_uniprot_mapping.items():
        target_aligned_sequence = seq_aligned[target_id]
        query_aligned_sequence = seq_aligned[query_id]
        target_aa_position_array = alignment_string_to_numbers(target_aligned_sequence)
        query_aa_position_array = alignment_string_to_numbers(query_aligned_sequence)
        all_mcs_positions = list(range(len(query_aa_position_array)))
        global_sequence_mappings[query_id] = dict(zip(query_aa_position_array, all_mcs_positions))
        global_sequence_mappings_reverse[query_id] = dict(zip(all_mcs_positions, query_aa_position_array))

        meta_uniprot_data["query_id"].append(query_id)
        meta_uniprot_data["target_id"].append(target_id)
        meta_uniprot_data["query_seq_aligned"].append(query_aligned_sequence)
        meta_uniprot_data["target_seq_aligned"].append(target_aligned_sequence)
        meta_uniprot_data["mapping"].append(dict(zip(query_aa_position_array, all_mcs_positions)))
        meta_uniprot_data["reverse_mapping"].append(dict(zip(all_mcs_positions, query_aa_position_array)))

    return pd.DataFrame(meta_uniprot_data)


def save_as_fasta(meta_uniprot_data_df, prefix = ""):
    for target_id in meta_uniprot_data_df.target_id.unique():
        fname = prefix + "alignment_uniprots_{}.fasta".format(target_id)
        subset_data = meta_uniprot_data_df[meta_uniprot_data_df.target_id == target_id]
        with open(fname, "w") as file:
            for row in subset_data.iterrows():
                target_seq = row[1].target_seq_aligned
                query_seq = row[1].query_seq_aligned
                query_id = row[1].query_id
                file.write(">"+target_id+"\n")
                file.write(target_seq+"\n")
                file.write(">"+query_id+"\n")
                file.write(query_seq+"\n")
        
        CrystalUtils(pdb_codes = [], dst_folder = ".").visualizeMSA(fname)

        

def get_topology_sequences(results_df):
    # load all structures into biopython tolopgy
    pdbparser = Bio.PDB.PDBParser(QUIET=True)   # suppress PDBConstructionWarning
    bio_pdb_topologies = []
    for row in results_df.iterrows():
        #print(row[1]["pdb_id"])
        pdb_id = row[1]["pdb_id"]
        struct = pdbparser.get_structure(pdb_id, row[1]['file_loc'])
        bio_pdb_topologies.append(struct)
    
    # extract per chain sequences
    pdb_file_seqs = {"pdb_id":[], 
                     "entity_id":[], 
                     "chain":[],
                     "file_seq":[],
                     "file_resnames":[],
                     "file_seq_ids":[]
                    }

    for i,topology in enumerate(bio_pdb_topologies):

        #details about the PDB
        pdb_id_details = results_df.iloc[i]
        pdb_id = pdb_id_details["pdb_id"]
        asym_id = pdb_id_details["first_asym_id"]
        entity_id = pdb_id_details["entity_id"]

        entity_details = entity_instance_mapping_df[(entity_instance_mapping_df["asym_id"]==asym_id) &
                                                (entity_instance_mapping_df["entity_id"]==entity_id)
                                                ]
        chain_id = entity_details.iloc[0].auth_asym_id

        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        print("PDB_ID: {}".format(pdb_id))
        print("ENTITY_ID: {}".format(entity_id))
        chains = topology.get_chains()
        for chain in chains:
            if chain.id == chain_id:
                print("chain: {}".format(chain.id))
                residues = list(chain.get_residues())
                res_names = [r.resname for r in residues]
                res_ids = [r.full_id[-1][-2] for r in residues]
                pdb_file_seqs["pdb_id"].append(pdb_id)
                pdb_file_seqs["entity_id"].append(entity_id)
                pdb_file_seqs["chain"].append(chain.id)
                pdb_file_seqs["file_seq"].append("".join([CrystalUtils.d3to1[r] for r in res_names]))
                pdb_file_seqs["file_resnames"].append(res_names)
                pdb_file_seqs["file_seq_ids"].append(res_ids)     
                
    return  pd.DataFrame(pdb_file_seqs)


def find_gaps(seq):
    return np.argwhere(np.array(list(seq))=="-").flatten()