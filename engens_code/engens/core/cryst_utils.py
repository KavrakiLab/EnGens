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
import os
from tabulate import tabulate 
import math
import numpy as np

# get PDBRenum files
def get_pdb_bioassembly_pdbrenum(pdb_id, dst):
    
    res_file = "{}/{}_renum.pdb1".format(dst, pdb_id)
    if not os.path.exists(res_file):
        cmd_get_pdbrenum = "wget http://dunbrack3.fccc.edu/PDBrenum/output_PDB_assembly/{}_renum.pdb1.gz".format(pdb_id)
        ret = os.system(cmd_get_pdbrenum)
        if ret != 0:
            return None
        cmd_get_pdbrenum = "mv {}_renum.pdb1.gz {}/{}_renum.pdb1.gz".format(pdb_id, dst, pdb_id)
        ret = os.system(cmd_get_pdbrenum)
        cmd_ungz = "gzip -d {}/{}_renum.pdb1.gz".format(dst, pdb_id)
        ret = os.system(cmd_ungz)
    else:
        print("Found existing {}".format(res_file))
    return res_file

# fix files with PDBFixer
def fix_pdb_files(input_file, output_file):
    if not os.path.exists(output_file):
        pdbfixer_command = "pdbfix.py -i {} -o {}".format(input_file,output_file)
        os.system(pdbfixer_command)
    else: 
        print("Found existing {}".format(output_file))
    return output_file

# for the given PDB entry get associated entities
def rscb_entities_from_entries(pdb_ids):
    pdb_ids_str = "[\""+"\",\"".join(pdb_ids)+"\"]"

    # QUERY - get the entity ids for the given entry ids
    query = """query {
      entries(entry_ids: """+pdb_ids_str+"""){
        polymer_entities {
          rcsb_id
          rcsb_polymer_entity_container_identifiers {
            reference_sequence_identifiers {
              database_accession
              database_name
            }
          }
        }
      }
    }"""  
    
    url = 'https://data.rcsb.org/graphql'
    r = requests.post(url, json={'query': query})
    print(r.status_code)
    
    if r.status_code == 200:
        json_data = r.json()

        #reformat dictionary
        ref_json_data = []
        for i, elem in enumerate(json_data["data"]["entries"]):
            pdb_id = pdb_ids[i]
            elem_polymers = elem["polymer_entities"]
            for ep in elem_polymers:
                rcsb_id = ep["rcsb_id"]
                if ep["rcsb_polymer_entity_container_identifiers"]["reference_sequence_identifiers"] is not None:
                    for seq_id in ep["rcsb_polymer_entity_container_identifiers"]["reference_sequence_identifiers"]:
                        assecion = seq_id["database_accession"]
                        db_name = seq_id["database_name"]
                        ref_json_data.append({"pdb_id":pdb_id, \
                                              "entity_id": rcsb_id,\
                                              "accession": assecion,\
                                              "database": db_name})
                else:
                    ref_json_data.append({"pdb_id":pdb_id, \
                                          "entity_id": rcsb_id,\
                                          "accession": None,\
                                          "database": None})
                    
        results_df = pd.DataFrame(ref_json_data)
        return results_df
    else:
        print("HTTP Request error {} for RSCB polymer_entities query".format(r.status_code))
        return None
    
    
#For the given entity ID - get the chain number
def rscb_polymer_chains_info(entity_ids):
    
    entity_id_str = "[\"" + "\",\"".join(entity_ids)+"\"]"
    
    # QUERY - get the entry ids for all the above to collect asym_ids
    query = """
    query {
      polymer_entities(entity_ids:"""+entity_id_str+""") {
        rcsb_id
        rcsb_entity_source_organism {
          ncbi_taxonomy_id
          ncbi_scientific_name
        }
        rcsb_cluster_membership {
          cluster_id
          identity
        }
        rcsb_polymer_entity_container_identifiers{
          asym_ids
        }
      }
    }
    """
    
    url = 'https://data.rcsb.org/graphql'
    r = requests.post(url, json={'query': query})
    print(r.status_code)
    
    if r.status_code == 200:
        json_data = r.json()
        polymer_entity_data = json_data['data']["polymer_entities"]
        entity_instance_connection = {"entity_id": [], "asym_ids":[]}
        for elem in polymer_entity_data:
            rcsb_id = elem["rcsb_id"]
            identifiers = elem["rcsb_polymer_entity_container_identifiers"]
            entity_instance_connection["asym_ids"].append(identifiers["asym_ids"])
            entity_instance_connection["entity_id"].append(rcsb_id)
        entity_instance_connection_df = pd.DataFrame(entity_instance_connection)
        return entity_instance_connection_df
    else:
        print("HTTP Request error {} for RSCB polymer_entities query".format(r.status_code))
        return None


def rscb_get_author_instance_info(instance_ids):
    
    entity_instance_ids = "[\"" + "\",\"".join(instance_ids)+"\"]"
    
    # QUERY - find the author chain id and author sequence mapping
    query = """
    query {
      polymer_entity_instances(instance_ids: """+entity_instance_ids+""") {
        rcsb_id
        rcsb_polymer_entity_instance_container_identifiers {
          asym_id
          auth_asym_id,
          auth_to_entity_poly_seq_mapping
        }
      }
    }

    """
    
    url = 'https://data.rcsb.org/graphql'
    r = requests.post(url, json={'query': query})
    print(r.status_code)
    
    if r.status_code == 200:
        json_data = r.json()
        polymer_entity_instances_data = json_data['data']["polymer_entity_instances"]
        entity_instances_mappings = {"instance_id": [], "asym_id":[], "auth_asym_id":[], "auth_seq_map":[]}
        for elem in polymer_entity_instances_data:
            rcsb_id = elem["rcsb_id"]
            identifiers = elem["rcsb_polymer_entity_instance_container_identifiers"]
            entity_instances_mappings["instance_id"].append(rcsb_id)
            entity_instances_mappings["asym_id"].append(identifiers["asym_id"])
            entity_instances_mappings["auth_asym_id"].append(identifiers["auth_asym_id"])
            entity_instances_mappings["auth_seq_map"].append(identifiers["auth_to_entity_poly_seq_mapping"])
        entity_instance_mapping_df = pd.DataFrame(entity_instances_mappings)
        return entity_instance_mapping_df
    else:
        print("HTTP Request error {} for RSCB polymer_entity_instances query".format(r.status_code))
        return None

def uniprot_get_details(uniprot_ids, include_pdb_info = False):
    
    uniprot_details = {"accession_id":[], 
                       "id":[],
                       "full_name":[],
                       "seq" : []}
    pdb_info = None
    
    if include_pdb_info:
        pdb_info = {"accession_id":[], 
                       "pdb_id":[],
                       "method":[],
                       "chains" : [],
                       "resolution" : []}
    uniprot_accession_url = "https://www.ebi.ac.uk/proteins/api/proteins/"
    
    for uni_id in uniprot_ids:
        accession_query = uniprot_accession_url+uni_id
        result_uniprot_details = requests.get(accession_query)
        if result_uniprot_details.status_code == 200:
            res_json = result_uniprot_details.json()
            uniprot_details["accession_id"].append( res_json["accession"] )
            uniprot_details["id"].append( res_json["id"] )
            if "protein" in res_json and "recommendedName" in res_json["protein"]:
                uniprot_details["full_name"].append( res_json['protein']['recommendedName']['fullName']['value'] )
            else:
                uniprot_details["full_name"].append("-")
              
            uniprot_details["seq"].append( res_json['sequence']['sequence'] )
            
            if include_pdb_info:
                for elem in res_json['dbReferences']:
                    if elem["type"] == "PDB":
                        pdb_info["accession_id"].append( res_json["accession"] )
                        pdb_info["pdb_id"].append( elem["id"] )
                        pdb_info["method"].append( elem["properties"]["method"] )
                        pdb_info["chains"].append( elem["properties"]["chains"] )
                        if "resolution" in elem["properties"]:
                            pdb_info["resolution"].append( elem["properties"]["resolution"] )
                        else:
                            pdb_info["resolution"].append("-")
        else:
            print("Uniprot query failed: response "+result_uniprot_details.status_code)
            return None
        
    return pd.DataFrame(uniprot_details), pd.DataFrame(pdb_info)


def get_pdb_metadata(pdb_ids):
    #pdb_id to entity (splitting chains, geting uniprot data)
    print("Fetching metadata for PDB entries - pdb id and related entity ids")
    entity_df = rscb_entities_from_entries(pdb_ids)
    #entity to instance (mapping chains to uniprot ids)
    print("Fetching metadata for PDB entries - entity ids and related uniprot ids")
    entity_ids = list(entity_df["entity_id"].unique())
    entity_instance_connection_df = rscb_polymer_chains_info(entity_ids)
    entity_instance_connection_df["first_asym_id"] = entity_instance_connection_df["asym_ids"].apply(lambda x: x[0])
    entity_df = entity_df.merge(entity_instance_connection_df, how="left", on="entity_id")
    entity_df["instance_id"] = entity_df["pdb_id"].apply(lambda x: x.upper())+"."+entity_df["first_asym_id"]
    instance_ids = list(entity_df["instance_id"].unique())
    #get author sequence mappings for each instance
    print("Fetching metadata for PDB entries - instance ids and related sequences")
    entity_instance_mapping_df = rscb_get_author_instance_info(instance_ids)
    entity_instance_mapping_df["entity_id"] = entity_instance_mapping_df.merge(entity_df, how="left", on="instance_id").entity_id
    entity_df = entity_df[~entity_df.accession.isna()]
    #get uniprot metadata
    print("Fetching metadata for UniProt ids")
    uniprot_ids = entity_df["accession"].unique()
    uniprot_details, _ = uniprot_get_details(uniprot_ids)
    return (entity_df, uniprot_details)

def map_uniprot_ids(uniprot_df):
    uniprot_ids = uniprot_df["accession_id"].unique()
    print("User mapping of the following UniProt ids to equivalent groups") 
    print(tabulate(uniprot_df[["accession_id", "id", "full_name"]]))
    possible_values = uniprot_ids
    final_mapping = {}
    print("Please map all entries to one of {}".format(possible_values))
    print("Or just press enter to map to self")
    for key in possible_values:
        key_row = uniprot_df[uniprot_df["accession_id"]==key]
        key_desc = key_row["full_name"]
        user_map = input("Map {} ({}) to:".format(key, key_desc.iloc[0]))
        if user_map is None or len(user_map)==0:
            print(key)
            user_map = key
        else:
            while not user_map in possible_values:
                print("Invalid input - try again")
                print("With one of {}".format(possible_values))
                user_map = input("Map {} ({}) to:".format(key, key_desc))
                if user_map is None or len(user_map)==0:
                    print(key)
                    user_map = key
        final_mapping[key] = user_map
        map_row = uniprot_df[uniprot_df["accession_id"]==user_map]
        map_desc = key_row["full_name"]
        print("Success!")
    print("Finall mapping:")
    print(final_mapping)
    return final_mapping

def input_domain_n(uniprot_n):
    while True: 
        domain_inp = input("Select the number of domains/chains in your complex: ")
        try:
            domain_n = int(domain_inp)
            if domain_n < 0:
                print("Number of domains/chains must be positive, try again")
            else:
                if domain_n > uniprot_n:
                    print("Number of domains/chains must be less than number of uniprot accessions in the dataset, try again")
                else:
                    return domain_n
        except:
            print("Number of domains/chains must be an integer, try again")
            
            
def input_domain_uniprot(domain_name, uniprot_ids):
    
    while True: 
        domain_inp = input("Input one of {} uniprot acession to map to {}: ".format(uniprot_ids, domain_name))
        try:
            domain_n = str(domain_inp)
            if not domain_n in uniprot_ids:
                print("You must select one of available uniprot accession ids: {}".format(uniprot_ids))
            else:
                return domain_n
        except:
            print("Number of domains/chains must be a string")

def pick_PDBs(pdb_df):
    print("Choose from: ")
    print(tabulate(pdb_df, \
                   headers = pdb_df.columns))
    print()
    while True:
        in_sele = input("Input indices of structures you want to include separated by space (e.g., '0 3 5' ) \n to select all input string all\n")
        try: 
            if "all" in in_sele:
                in_list = list(range(pdb_df.shape[0]))
            else:
                in_list = [int(x) for x in in_sele.split(" ")]
        except:
            print("Invalid input, try again")
            continue

        print(in_list)
        correct_selection = True
        for elem in in_list:
            if not elem in list(range(pdb_df.shape[0])):
                correct_selection = False
                print("Incorrect input, try again")
                print(pdb_df.shape)
                break
        if not correct_selection:
            continue
        selected_pdbs = pdb_df.iloc[in_list]
        print("Selected PDBs:")
        print(tabulate(selected_pdbs, \
                       headers = selected_pdbs.columns))
        return selected_pdbs
        break
    
    
def define_domains(uniprot_df, entity_df):
    uniprot_ids = uniprot_df["accession_id"].unique()
    
    # define the number of domains
    domain_n = input_domain_n(len(uniprot_ids))
    print("Total of #{} domains!".format(domain_n))
    print()
    
    # map each accession ID to a domain
    print("Map each domain to uniprot accession (by their accession_id).")
    uniprot_domain_map = {}
    domain_uniprot_map = {"domain_name":[], "uniprot_id":[]}
    selected_uniprots = []
    for i in range(domain_n):
        available_uniprots = uniprot_df[~uniprot_df.accession_id.isin(selected_uniprots)].reset_index()
        domain_name = "DOMAIN{}".format(i)
        print()
        print("Map uniprot accessions to {}".format(domain_name))
        print("Choose from: ")
        print(tabulate(available_uniprots[["accession_id", "id", "full_name"]], \
                       headers = ["accession_id", "id", "full_name"]))
        print()
        while True:
            in_sele = input("Input indices of uniprot metadata separated by space (e.g., '0 3 5' ) \n to select all input string all\n")
            try: 
                if "all" in in_sele:
                    in_list = list(range(available_uniprots.shape[0]))
                else:
                    in_list = [int(x) for x in in_sele.split(" ")]
            except:
                print("Invalid input, try again")
                continue
                
            print(in_list)
            correct_selection = True
            for elem in in_list:
                if not elem in list(range(available_uniprots.shape[0])):
                    correct_selection = False
                    print("Incorrect input, try again")
                    print(available_uniprots.shape)
                    break
            if not correct_selection:
                continue
            selected_unis = available_uniprots.iloc[in_list]
            print("Selected for DOMAIN{}: ".format(i))
            print(tabulate(selected_unis[["accession_id", "id", "full_name"]], \
                           headers = ["accession_id", "id", "full_name"]))
            selected_acc_ids = list(selected_unis.accession_id)
            selected_uniprots.extend(selected_acc_ids)
            for uni_acc_id in selected_acc_ids:
                uniprot_domain_map[uni_acc_id] = domain_name
                domain_uniprot_map["domain_name"].append(domain_name)
                domain_uniprot_map["uniprot_id"].append(uni_acc_id)
            break
    print()
    #check if all accession IDs were mapped!
    for uniprot_id in uniprot_ids:
        if not uniprot_id in uniprot_domain_map:
            print("WARNING: uniprot id {} not mapped to any domain".format(uniprot_id))
            uniprot_domain_map[uniprot_id] = None
            
    uniprot_df["domain"] = uniprot_df["accession_id"].apply(lambda x: uniprot_domain_map[x])
    entity_df["domain"]=entity_df["accession"].apply(lambda x: uniprot_domain_map[x])
    domain_df = pd.DataFrame(domain_uniprot_map)
    return uniprot_df, entity_df, domain_df

from Bio.PDB import Select, PDBIO
from Bio.PDB.PDBParser import PDBParser

class ChainSelect(Select):
    def __init__(self, chain):
        self.chain = chain

    def accept_chain(self, chain):
        if chain.get_id() == self.chain:
            return 1
        else:          
            return 0

def extract_chain(in_file, chain, out_file):
    p = PDBParser(PERMISSIVE=1, QUIET=True)       
    structure = p.get_structure(in_file, in_file)
    pdb_chain_file = out_file    
    io_w_no_h = PDBIO()               
    io_w_no_h.set_structure(structure)
    io_w_no_h.save('{}'.format(pdb_chain_file), ChainSelect(chain))

def extract_chains_create_file(file_location, domain, accession, chain):
    new_file_dir = Path(path.join(Path(file_location).parent, domain+"/"+accession))
    new_file_dir.mkdir(parents=True, exist_ok=True)
    new_file = path.join(new_file_dir, Path(file_location).stem+"_"+chain+".pdb")
    print(new_file)
    extract_chain(file_location, chain, new_file)
    try:
        tmp = pd.read_csv(new_file, sep="\n", header=None)
        if tmp.shape[0] < 3:
            raise None
    except:
        print("Warning: File {} does not contain chain {}".format(file_location, chain))
        os.remove(new_file)
        return None
    return new_file

def pdb_reorder_residues(file):
    print(file)
    pd.options.mode.chained_assignment = None
    tmp = pd.read_csv(file, sep="\n", header=None)
    tmp_atom = tmp[tmp[0].apply(lambda x: x.split(" ")[0]=="ATOM")]
    tmp_atom["len"] = tmp_atom[0].apply(lambda x: len(x.split()))
    tmp_atom["res_n"] = tmp_atom[0].apply(get_residue_number)
    tmp_atom["res_id_3"] = tmp_atom[0].apply(get_residue_id)
    tmp_atom["res_id_1"] = tmp_atom["res_id_3"].apply(lambda x: d3to1[x])
    tmp_atom.sort_values(by=["res_n"], inplace=True)
    tmp_atom[[0]].to_csv(file, sep="\n", index=False, header=False)
    
d3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
        'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
        'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
        'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M', 'YCM':'C'}

def get_residue_number(line):
    split_line = line.split()
    if len(split_line) == 11 and len(split_line[4])>1:
        return int(split_line[4][1:])
    elif len(split_line) == 12:
        return int(split_line[5])
    else:
        return int(split_line[5])
    
def get_residue_id(line):
    split_line = line.split()
    return split_line[3]


def get_max_residue_n(file):
    pd.options.mode.chained_assignment = None
    tmp = pd.read_csv(file, sep="\n", header=None)
    #print("x")
    tmp_atom = tmp[tmp[0].apply(lambda x: x.split(" ")[0]=="ATOM")]
    #print("x")
    tmp_atom["len"] = tmp_atom[0].apply(lambda x: len(x.split()))
    #print("x")
    tmp_atom["res_n"] = tmp_atom[0].apply(get_residue_number)
    return tmp_atom["res_n"].max()

def get_min_residue_n(file):
    pd.options.mode.chained_assignment = None
    tmp = pd.read_csv(file, sep="\n", header=None)
    #print("x")
    tmp_atom = tmp[tmp[0].apply(lambda x: x.split(" ")[0]=="ATOM")]
    #print("x")
    tmp_atom["len"] = tmp_atom[0].apply(lambda x: len(x.split()))
    #print("x")
    tmp_atom["res_n"] = tmp_atom[0].apply(get_residue_number)
    return tmp_atom["res_n"].min()

def extract_msa_sequence(file, min_len, max_len):
    
    tmp = pd.read_csv(file, sep="\n", header=None)
    tmp_atom = tmp[tmp[0].apply(lambda x: x.split(" ")[0]=="ATOM")]
    #tmp_atom["len"] = tmp_atom[0].apply(lambda x: len(x.split()))
    tmp_atom["res_n"] = tmp_atom[0].apply(get_residue_number)
    tmp_atom["res_id_3"] = tmp_atom[0].apply(get_residue_id)
    tmp_atom["res_id_1"] = tmp_atom["res_id_3"].apply(lambda x: d3to1[x])
    tmp_atom = tmp_atom.sort_values(by=["res_n"])
    tmp_atom = tmp_atom.groupby(by=["res_n", "res_id_1"]).count().reset_index()
    mask = pd.DataFrame({"res_n": list(range(min_len, max_len+1)),
                        "res_id_1": ["-" for i in range(min_len, max_len+1)]})
    marged_df = mask.merge(tmp_atom, on="res_n", how='left')
    full_string = list(marged_df.apply(lambda x: x.res_id_1_y,axis = 1))
    full_string = ["-" if (type(x)==float and math.isnan(x)) else str(x) for x in full_string ]
    #print(full_string)
    full_string = "".join(full_string)
    return full_string

# Extract MSA for the domain 
# Based on the residue numbers (for the elements that have same uniprot id source)
def extract_msa_for_domain_signle_uniprot(domain, uniprot_id, msa_result_file, entity_df):
    
    print("Creating MSA for domain {} based on {} uniprot accession".format(domain, uniprot_id))

    tmp_df = entity_df[(entity_df.domain==domain) & (entity_df.accession==uniprot_id)]
    #print(tabulate(tmp_df, \
    #               headers = tmp_df.columns))
    start_res = tmp_df.domain_file_loc.apply(lambda x: get_min_residue_n(x)).min()
    end_res = tmp_df.domain_file_loc.apply(lambda x: get_max_residue_n(x)).max()
    tmp_df["seq"] = tmp_df.domain_file_loc.apply(lambda x: extract_msa_sequence(x, start_res, end_res))
    msa = [">{}\n{}\n".format(row.instance_id, row.seq) for i, row in tmp_df.iterrows()]
    msa = "".join(msa)
    with open(msa_result_file, "w") as in_file:
        in_file.write(msa)
    print()
    print("Extracted MSA for domain {} based on {} uniprot accession to {}".format(domain, uniprot_id, msa_result_file))
    print("Residue range {}-{}".format(start_res, end_res))
    print("Check out the MSA at {}".format(msa_result_file))
    return start_res, end_res, msa_result_file
    #files = list(tmp_df.domain_file_loc.unique())

def extract_mcs(pdb_files, mcs, start, out_dir = "."):
    new_files = []
    for file in os.listdir(pdb_files):
        tmp = pd.read_csv(path.join(pdb_files, file), sep="\n", header=None)
        tmp["len"] = tmp[0].apply(lambda x: len(x.split()))
        tmp["res_n"] = tmp[0].apply(get_residue_number)
        tmp["res_id_3"] = tmp[0].apply(get_residue_id)
        tmp["res_id_1"] = tmp["res_id_3"].apply(lambda x: d3to1[x])
        tmp = tmp[tmp.res_n.isin(mcs+start)].reset_index()
        new_file = path.join(out_dir, path.join(pdb_files, Path(file).stem+"_mcs.pdb"))
        tmp[0].to_csv(new_file, sep="\n", index=False, header=False)
        new_files.append(new_file)
    return new_files

def check_atom_numbers_mcs_match(pdb_files):
    all_mcs_df = None
    for file in os.listdir(pdb_files):
        if Path(file).stem[-3:] == "mcs":
            tmp = pd.read_csv(path.join(pdb_files, file), sep="\n", header=None)
            tmp["res_n"] = tmp[0].apply(get_residue_number)
            tmp["res_id_3"] = tmp[0].apply(get_residue_id)
            tmp["atom"] = tmp[0].apply(lambda x: x.split()[2])
            tmp["res_id_1"] = tmp["res_id_3"].apply(lambda x: d3to1[x])
            tmp["file"] = [file for i in range(tmp.shape[0])]
            if all_mcs_df is None:
                all_mcs_df = tmp
            else:
                all_mcs_df = pd.concat([all_mcs_df, tmp])
    all_mcs_df_backbone = all_mcs_df[all_mcs_df.atom.isin(["N", "CA", "C", "O"])]
    all_mcs_df_backbone.groupby(by=["file"]).count().reset_index()
    if len(list(all_mcs_df_backbone[0].unique())) == 1:
        print("ERROR: backbones of MCS for domain {} asseccion {} do not match - discarding this domain")
        return False
    return True

import itertools

def intervals_extract(iterable):
    iterable = sorted(set(iterable))
    for key, group in itertools.groupby(enumerate(iterable),
    lambda t: t[1] - t[0]):
        group = list(group)
        yield [group[0][1], group[-1][1]]

def get_mcs_residues_from_fasta(fasta_file, start_res=0):
    ali_df = pd.read_csv(fasta_file, comment=">", header=None, sep="\n")
    ali_df = ali_df[0].apply(lambda x: [False if x== "-" else True for x in list(x)])
    ali_np = np.array(ali_df.to_list())
    ali_np_all = np.all(ali_np, axis=0)
    mcs = np.argwhere(ali_np_all).flatten()+start_res
    return(mcs)

def load_fasta(fasta):
    tmp = pd.read_csv(fasta, sep="\n", header=None)
    names = []
    seqs = []
    curr_seq = []
    for i, line in tmp.iterrows():
        line=line[0]
        if line[0]==">":
            names.append(line)
            if len(curr_seq)>0:
                seqs.append("".join(curr_seq))
                curr_seq = []
        else:
            curr_seq.append(line)
    seqs.append("".join(curr_seq))
    return pd.DataFrame({"name": names, "seq": seqs})

def get_resn_from_pdb(orig_file):
    tmp = pd.read_csv(orig_file, sep="\n", header=None)
    tmp_atom = tmp[tmp[0].apply(lambda x: x.split(" ")[0]=="ATOM")]
    tmp_atom["len"] = tmp_atom[0].apply(lambda x: len(x.split()))
    tmp_atom["res_n"] = tmp_atom[0].apply(get_residue_number)
    return list(tmp_atom.groupby(by="res_n").count().reset_index().res_n)

def map_mcs_to_resn(seq_msa, mcs_all, orig_file):
    #instance: 1-N for each caracter in MSA (N size of MSA)
    #instance_indices: -1 for gap 0-N for other residues
    i=0
    instance_indices = []
    for elem in seq_msa:
        if elem =="-":
            instance_indices.append(-1)
        else:   
            instance_indices.append(i)
            i+=1
    instance_indices = np.array(instance_indices)
    mcs_instance = instance_indices[mcs_all]
    #local: 1-n for each letter in MSA (n size of this proten)
    local_indices = list(range(i))
    #original: xxx for resn in the files - correspond to new local indices
    resn_indices = get_resn_from_pdb(orig_file)

    resn_local_mapping = dict(zip(resn_indices, local_indices))
    local_resn_mapping = dict(zip(local_indices, resn_indices))

    mcs_instance_resn = [local_resn_mapping[elem] for elem in mcs_instance]
    return mcs_instance_resn


def perform_mTM_align(domain_metadata, pdb_metadata, dst_structure):
    
    grouped_domains = domain_metadata.groupby(by=["domain_name"]).count().reset_index()
    multi_ac_domains = grouped_domains[grouped_domains.mcs>1].domain_name.unique()
    print("Multi-accession domains (to be processed): {}".format(multi_ac_domains))

    multi_uni_mcs_metadata = {"domain":[], "uniprot_id":[], \
                              "template_file": [], "mcs_resn":[], \
                              "mcs_resn_regions":[]}
    
    if len(multi_ac_domains) == 0:
        return None
    
    for domain in multi_ac_domains:
        print("Processing domain {}".format(domain))
        root_dir = Path(path.join(dst_structure, domain+"/mTM_align"))
        root_dir.mkdir(parents=True, exist_ok=True)
        # Create mTM-align input 
        mtm_input_file = path.join(root_dir, "mTM_input.txt")
        tmp_df = domain_metadata[domain_metadata.domain_name==domain]
        uniprot_ids = list(tmp_df.uniprot_id)

        #select an single mcs structure from each dom-uni pair for mTM-align
        files = []
        mcs_files = []
        sample_pdbs = []
        for uniprot_id in uniprot_ids:
            sample_pdb = pdb_metadata[(pdb_metadata.accession == uniprot_id) &\
                                               (pdb_metadata.domain == domain)].iloc[0]
            mcs_file = path.join(Path(sample_pdb.domain_file_loc).parent, Path(sample_pdb.domain_file_loc).stem+"_mcs.pdb")

            mcs_files.append(mcs_file)
            sample_pdbs.append(sample_pdb.pdb_id)
            cp_loc_file = path.join(root_dir,Path(sample_pdb.domain_file_loc).stem+"_mcs.pdb")
            os.system("cp {} {}".format(mcs_file, cp_loc_file))
            files.append(cp_loc_file)
        
        mtm_al_metadata = pd.DataFrame({"file": files, "mcs_file": mcs_files, \
                           "sample_pdb_id": sample_pdbs, "uniprot_id": uniprot_ids,\
                          "domain": [domain for i in range(len(files))]})
        
        with open(mtm_input_file, "w") as inf:
            inf.write("\n".join(files))

        #run mTM-align
        print()
        print("Running mTM-align of {} for uniprots {} ".format(domain, uniprot_ids))
        system("mTM-align -i "+mtm_input_file
           +" -outdir {}".format(root_dir) 
           +" > mTM-progress.txt")
        
        result_mTM_fasta = path.join(root_dir, "result.fasta")
        #generate mTM-alignment for all files
        fasta_df = load_fasta(result_mTM_fasta)
        fasta_2d = np.array(fasta_df.seq.apply(lambda x: list(x)).to_list())
        bool_2d = np.array(fasta_df.seq.apply(lambda x: [False if c=="-" else True for c in list(x)]).to_list())
        mcs_all = np.argwhere(bool_2d.all(axis=0)).flatten()

        #map MCS to each instance
        for i, row in mtm_al_metadata.iterrows():
            fasta_title = ">"+str(Path(row.file).stem)+".pdb"
            file = row.mcs_file
            file_msa_seq = fasta_df[fasta_df["name"] == fasta_title].seq.iloc[0]
            mcs_resn = map_mcs_to_resn(file_msa_seq, mcs_all, file)
            mcs_resn_regions = [tmp for tmp in intervals_extract(mcs_resn)]
            multi_uni_mcs_metadata["domain"].append(domain)
            multi_uni_mcs_metadata["uniprot_id"].append(uniprot_ids[i])
            multi_uni_mcs_metadata["template_file"].append(file)
            multi_uni_mcs_metadata["mcs_resn"].append(mcs_resn)  
            multi_uni_mcs_metadata["mcs_resn_regions"].append(mcs_resn_regions)  
    return pd.DataFrame(multi_uni_mcs_metadata)

def find_pdb_ids_dom_uni_qualifying(pdb_metadata, selected_unis):
    
    selected_pdb_uniprot_ids = {"pdb_id": [], "domain":[], "uniprot_id":[]}
    for pdb_id in list(pdb_metadata.pdb_id.unique()):
        tmp_df = pdb_metadata[pdb_metadata["pdb_id"]==pdb_id]
        #needs to have all selected domains from the selected accession 
        selected_domains = list(selected_unis.domain_name.unique())
        qualifies = True

        pdb_ids = []
        domains = []
        accessions = []
        for dom in selected_domains:
            related_unis = list(selected_unis[selected_unis.domain_name == dom].uniprot_id.unique())
            selected_df = tmp_df[(tmp_df.domain==dom) & (tmp_df.accession.isin(related_unis))]
            if len(selected_df) == 0:
                qualifies = False
                break
            else:
                pdb_ids.append(pdb_id)
                domains.append(dom)
                accessions.append(selected_df.accession.iloc[0])
        if qualifies:
            selected_pdb_uniprot_ids["pdb_id"].extend(pdb_ids)
            selected_pdb_uniprot_ids["domain"].extend(domains)
            selected_pdb_uniprot_ids["uniprot_id"].extend(accessions)

    return pd.DataFrame(selected_pdb_uniprot_ids)

def visualizeMSA(fasta_file, regions=[], template_loc=None):

    # Read in the resulting MSA
    mulseq_algn = []
    names = []
    curr_line = None
    fasta_seqs = []
    
    tmp_df = pd.read_csv(fasta_file, header=None, sep="\n")
    mulseq_algn = list(tmp_df[tmp_df[0].apply(lambda x: x[0] is not ">")][0])
    names = list(tmp_df[tmp_df[0].apply(lambda x: x[0] == ">")][0])
    for i, name in enumerate(names):
        fasta_seqs.append(FastaSequence(name, mulseq_algn[i]))
    
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
    if template_loc == None:
        #default template location
        template_loc = "{}/template_msa.html".format(str(Path(sys.executable).parent.parent))
                                          
    with open(template_loc, "rt") as fin:
        with open(output_name, "wt") as fout:
            for line in fin:
                new_line = line.replace('$VARIABLES_PLACEHOLDER$', js_txt)
                new_line = new_line.replace('$SL$', "\"{}\"".format(seq_length))
                fout.write(new_line)
    return output_name