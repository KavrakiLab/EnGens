import pandas as pd
import requests

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
        
def uniprot_get_details(uniprot_ids):
    uniprot_details = {"accession_id":[], 
                       "id":[],
                       "full_name":[],
                       "seq" : []}
    uniprot_accession_url = "https://www.ebi.ac.uk/proteins/api/proteins/"
    
    for uni_id in uniprot_ids:
        accession_query = uniprot_accession_url+uni_id
        result_uniprot_details = requests.get(accession_query)
        if result_uniprot_details.status_code == 200:
            res_json = result_uniprot_details.json()
            uniprot_details["accession_id"].append( res_json["accession"] )
            uniprot_details["id"].append( res_json["id"] )
            uniprot_details["full_name"].append( res_json['protein']['recommendedName']['fullName']['value'] )
            uniprot_details["seq"].append( res_json['sequence']['sequence'] )
            
        else:
            print("Uniprot query failed: response "+result_uniprot_details.status_code)
            return None
    return pd.DataFrame(uniprot_details)


def rscb_1d3d_get_residue_maps(results_df, meta_uniprot_mapping):
    
    query_part1 = """query {
      alignment(
        from:PDB_ENTITY,
        to:UNIPROT,
        queryId: " """
    query_part2 = """ " 
      ){
        query_sequence
        target_alignment {
          target_id
          target_sequence
          coverage{
            query_coverage
            query_length
            target_coverage
            target_length
          }
          aligned_regions {
            query_begin
            query_end
            target_begin
            target_end
          }
        }
      }
    }"""
    
    url = 'https://1d-coordinates.rcsb.org/graphql'

    mapping_1d_3d = {
        "entity_id":[],
        "uniprot_id":[],
        "seq_map":[],
        "seq_map_reverse":[],
        "seq":[],
        "target_seq":[],
        "coverage":[],
        "regions":[]
    }

    local_alignment_map = {}
    local_alignment_map_reverse = {}
    query_sequences = {}
    coverages = {}

    for i, row in results_df.iterrows():
        pdb_id = row["pdb_id"]
        entity_id = row["entity_id"]
        uniprot_id = row["accession"]
        meta_uniprot_id = meta_uniprot_mapping[uniprot_id]
        r = requests.post(url, json={'query': query_part1[:-1]+row["entity_id"]+query_part2[1:]})
        if r.status_code == 200:
            r_js = r.json()
            query_sequence = r_js["data"]["alignment"]["query_sequence"]
            N_alignments_per_entity = len(r_js["data"]["alignment"]["target_alignment"])
            for target_alignment in r_js["data"]["alignment"]["target_alignment"]:
                target_id = target_alignment["target_id"]
                target_cov = target_alignment["coverage"]["target_coverage"]
                target_len = target_alignment["coverage"]["target_length"]

                algn_regs = target_alignment["aligned_regions"]

                q_indices = []
                t_indices = []

                for region in algn_regs:
                    qb = int(region["query_begin"])-1
                    qe = int(region["query_end"])
                    tb = int(region["target_begin"])-1
                    te = int(region["target_end"])
                    q_list = list(range(qb, qe))
                    t_list = list(range(tb, te))
                    q_indices.extend(q_list)
                    t_indices.extend(t_list)

                chain = "-".join([pdb_id, entity_id, target_id])

                query_sequences[chain] = None
                if N_alignments_per_entity == 1:
                    query_sequences[chain] = (query_sequence, None, None)
                else:
                    query_sequences[chain] = (query_sequence, qb, qe)

                local_alignment_map[chain] = dict(zip(q_indices, t_indices))
                local_alignment_map_reverse[chain] = dict(zip(t_indices, q_indices))
                mapping_1d_3d["entity_id"].append(entity_id)
                mapping_1d_3d["uniprot_id"].append(target_id)
                mapping_1d_3d["seq_map"].append(dict(zip(q_indices, t_indices)))
                mapping_1d_3d["seq_map_reverse"].append(dict(zip(t_indices, q_indices)))
                mapping_1d_3d["seq"].append(query_sequence)
                mapping_1d_3d["target_seq"].append(target_alignment["target_sequence"])
                mapping_1d_3d["coverage"].append(target_cov)
                mapping_1d_3d["regions"].append(algn_regs)
    return pd.DataFrame(mapping_1d_3d)
