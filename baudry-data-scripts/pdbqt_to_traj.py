#!/usr/bin/env python
import subprocess
import os
import mdtraj
from pdbtools import pdb_delelem
from pdbtools import pdb_reatom

if __name__ == '__main__':
    import sys
    import getopt


    def usage():
        print ("Usage: pdbqt_to_traj.py -f foldername")
        print ("    Description of command...")
        print ("        [-f]    pdbqt_folder")
        print ("        [-s]    autodock prbqt_to_pdb.py script location")
        print ("        [-o]    trajname")
        print ("    Optional parameters:")
        print ("        [-v]    verbose output")


    # process command arguments
    try:
        opt_list, args = getopt.getopt(sys.argv[1:], 'f:s:o:v')

    except (getopt.GetoptError, msg):
        print ('pdbqt_to_traj.py: %s' %msg)
        usage()
        sys.exit(2)

    # initialize required parameters
    #-s: pdbqt_filename_stem
    pdbqt_foldername =  None
    pdbqt_script = None

    # optional parameters
    verbose = None
    trajname =  None 

    #parse input
    for o, a in opt_list:
        if o in ('-f', '--f'):
            pdbqt_foldername = a
            if verbose: print ('set pdbqt_foldername to ', pdbqt_foldername)
        if o in ('-s', '--s'):
            pdbqt_script = a
            if verbose: print ('set pdbqt_script to ', pdbqt_script)
        if o in ('-o', '--o'):
            trajname = a 
            if verbose: 
                print ('set output trajname to ', a)
        if o in ('-v', '--v'):
            verbose = True
            if verbose: print ('set verbose to ', True)
        if o in ('-h', '--'):
            usage()
            sys.exit()

    if not pdbqt_foldername:
        print ('pdbqt_to_pdb: pdbqt_foldername must be specified.')
        usage()
        sys.exit()  

    if not trajname:
        print ('trajname: trajname must be specified.')
        usage()
        sys.exit() 


    #open folder containing  
    files = []
    for filename in os.listdir(pdbqt_foldername):
        if filename.endswith(".pdbqt"): 
            #print(os.path.join(pdbqt_foldername, filename))
            files.append(os.path.join(pdbqt_foldername, filename))
            continue
        else:
            continue  

    if verbose: print(str(len(files)) + " pdbqt files to convert.")

    pdb_folder=pdbqt_foldername[:-1]+"_pdb"
    subprocess.run(["mkdir", "-p", pdb_folder])
    print(pdb_folder)
    pdb_files = []
    for file in files:
        pdb_name = pdb_folder+file[file.rfind("/"):file.find(".pdbqt")]+".pdb"
        pdb_noH_name = pdb_folder+file[file.rfind("/"):file.find(".pdbqt")]+"_noH.pdb"
        pdb_index = pdb_name[pdb_name.find("at-frame")+len("at-frame"):pdb_name.find(".pdb")]
        pdb_index = int(pdb_index)
        subprocess.call(["python2", pdbqt_script,  "-f", file, "-o", pdb_name, "-v"])
        pdb_noH = pdb_delelem.delete_elements(open(pdb_name), ["H"])
        with open("tmp.pdb", 'w') as f:
            for x in pdb_noH:
                f.write(str(x))
        pdb_noH = pdb_reatom.renumber_atom_serials(open('tmp.pdb'), 1)
        with open(pdb_noH_name, 'w') as f:
            for x in pdb_noH:
                f.write(str(x))
        subprocess.call(["rm", pdb_name])
        subprocess.call(["rm", "tmp.pdb"])
        
        pdb_files.append((pdb_noH_name, pdb_index))

    pdb_files.sort(key=lambda x: x[1])
    print(pdb_files)
    
    trajs = []
    for file, index in pdb_files:
        trajs.append(mdtraj.load(file))

    res_traj = None
    tmp_traj = []
    for i, traj in enumerate(trajs):
        if len(tmp_traj)>0:
            if not traj.n_atoms == tmp_traj[-1].n_atoms:
                print("Unequal atoms for trajectory: "+str(i*200)+" skipping")
                continue
        tmp_traj.append(traj)

    res_traj = mdtraj.join(tmp_traj)
    res_traj.save(trajname)

    
