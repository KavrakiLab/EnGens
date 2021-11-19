import os
import pandas as pd
import urllib.request
from Bio.PDB import PDBParser

#----------------Get sequence from pdb-----------------------
# You can use a dict to convert three letter code to one letter code
d3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M', 'YCM':'C'}

d1to3 = {v: k for k, v in d3to1.items()}

def getSequence(pdb_file):
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
            seq = []
            for residue in chain:
                if not residue.resname in d3to1:
                    break
                seq.append(residue.resname+str(residue.get_id()[1]))
            seqs.append(seq)
    return seqs
#--------------------------------

#================STEP 1 get sequence of the reference pdb======================================
folder_interactions = "./GCPCR_Ligands/GCPCR_Ligands/OPRK1/"

pdb_ref = "./oprk1_tested_pdb/at-frame0_noH.pdb"
ref_seq = getSequence(pdb_ref)[0]
#print(ref_seq)



#================STEP 2 read all listed interactions ======================================
ligand_interactions={}
for i, file in enumerate(os.listdir(folder_interactions)):
    #print(i)
    #find the pdb code
    pdb_code = file[file.rfind("_")+1: file.rfind("_")+5]
    #print(pdb_code)
    if not file[:len("Interaction")] =="Interaction": continue
    #get the info 
    info = pd.read_excel(folder_interactions+file, engine='openpyxl')
    info["aa-pos"] = info["Amino Acid"].apply(lambda x: d1to3[x]) + info["Sequence Number"].astype(str)
    info_aa_pos = info["aa-pos"]
    unique_aa_pos = pd.unique(info_aa_pos)
    #print(unique_aa_pos)
    ligand_interactions[pdb_code] = unique_aa_pos

#================STEP 3 fetch PDBs and get the sequences in the files======================================
#get sequences from PDBs
pdb_seqs = {}
for pdb_code, interations in ligand_interactions.items():
    pdb_download = folder_interactions+pdb_code+'.pdb'
    if not os.path.exists(pdb_download):
        urllib.request.urlretrieve('http://files.rcsb.org/download/'+pdb_code+'.pdb', pdb_download)
    seq = getSequence(pdb_download)
    pdb_seqs[pdb_code] = seq

#================STEP 4 find the sequence from each pdb that correctly corresponds to the interactions ======================================

def validate_seq(seq, interactions):
    for name_pos in interactions:
        if not name_pos in seq:
            return False

    return True

pdb_interaction_seq = {}
for pdb_code, seqs in pdb_seqs.items():
    interaction = ligand_interactions[pdb_code]
    found_seq = False
    for s in seqs:
        s_rev = [d3to1[x[:3]]+x[3:] for x in s]
        s = [d1to3[x[0]]+x[1:] for x in s_rev]
        if validate_seq(s, interaction):
            pdb_interaction_seq[pdb_code] = (interaction, s)
            found_seq = True
            break
    if not found_seq: 
        print("Didn't find the sequence with given interaction residues for:")
        print(pdb_code)
        print(interaction)
        print(seqs)

#print(pdb_interaction_seq)

#================STEP 5 align ref and pdb sequences  ==========================
#============== find the interaction residues in the reference structur ======================================
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Align import substitution_matrices
blosum62 = substitution_matrices.load("BLOSUM62")

def to_one_letter_seq(seq):
    res = ""
    for elem in seq:
        res+=d3to1[elem[:3]]
    return res



binding_site_residues = []

#rederence PDB to one-letter-code sequence
ref_seq_ol = to_one_letter_seq(ref_seq)
#print(ref_seq_ol)

#go through all ligand interaction pdb files
for pdb_code, elem in pdb_interaction_seq.items():
    #print("///////////////////////////////// PDB CODE //////////////////////////////////")
    #print(pdb_code)
    #print("\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\")
    interactions = elem[0]
    pdb_sequence = elem[1]
    pdb_sequence_ol = to_one_letter_seq(pdb_sequence)

    #do the alignment ref to pdb
    alignments = pairwise2.align.localds(ref_seq_ol, pdb_sequence_ol, blosum62, -10, -1)
    interaction_mapping=None

    min_failed = len(interactions)
    failed = None
    #find the alignment that maps the interactions correctly
    for al_num, alignment in enumerate(alignments):
        #print("######################ALIGNMENT "+str(al_num)+" #######################")
        #print(format_alignment(*alignment))
        #print("#########################################################################")
        seqA = alignment.seqA
        seqB = alignment.seqB
        ref_interaction_inds = []
        failed_interactions = []
        #check if all interactions are mapped
        for interaction in interactions:
            #print("---------------mapping interaction-----------------")
            #print(interaction)
            #print("--------------------------------------------------------")
            pdb_index = pdb_sequence.index(interaction)
            if pdb_index < 0: break
            seqB_ind = None
            cnt_no_gap = 0
            for cnt, char in enumerate(seqB):
                if char=="-": continue
                if cnt_no_gap == pdb_index:
                    seqB_ind = cnt
                    break
                cnt_no_gap+=1
            if seqB_ind == None:
                #print("unable to map interaction")
                break
            #check if there is the same letter at the positions
            if seqA[seqB_ind] == '-':
                #print("not mapped")
                #print((interaction, seqB_ind, seqB[seqB_ind], seqA[seqB_ind]))
                failed_interactions.append((interaction, seqB_ind, seqB[seqB_ind], seqA[seqB_ind], al_num))
                continue
            seqA_ind = seqB_ind
            ref_ind_ol = None
            cnt_no_gap = 0
            for cnt_gap, char in enumerate(seqA):
                if char == "-": continue
                if cnt_gap == seqA_ind:
                    ref_ind_ol = cnt_no_gap
                    break
                cnt_no_gap += 1
            #print("mapped")
            #print((interaction, seqB_ind, seqB[seqB_ind], seqA_ind, seqA[seqA_ind], ref_seq[ref_ind_ol]))
            ref_interaction_inds.append(ref_ind_ol)
        if len(failed_interactions) < min_failed:
            min_failed = len(failed_interactions)
            failed = failed_interactions
            interaction_mapping = (ref_interaction_inds, failed, al_num)
    
    print("///////////////////////////////// RESULTS FOR //////////////////////////////////")
    print(pdb_code)
    print("\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\")
    if interaction_mapping == None:
        print("COuldn't map any interactions any of the alignments to ref")
        print(pdb_code)
    else:
        print("Mapped interactions to alignment")
        print(format_alignment(*alignments[interaction_mapping[2]]))
        print("Failed to map:")
        print(str(len(interaction_mapping[1]))+"/"+str(len(interactions)))
        print("listing ------")
        for f in interaction_mapping[1]:
            print(f)
        print("------")
        binding_site_residues.append(interaction_mapping)
   
print(binding_site_residues)

residue_indices = []
residue_fullnames = []
for elem in binding_site_residues:
    residue_indices.extend(elem[0])
    for i in elem[0]:
        residue_fullnames.append(ref_seq[i])
print("--------------------------------------------------------------")
print(sorted(set(residue_indices)))
print(sorted(set(residue_fullnames), key=lambda x: int(x[3:])))