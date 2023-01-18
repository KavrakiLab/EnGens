
from openmm.app import Modeller
from pdbfixer import PDBFixer
from openmm.app import PDBFile

import argparse

def fix_file(input_file, out_file):
    fixer = PDBFixer(filename=input_file)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    print(fixer.missingAtoms)
    fixer.addMissingAtoms()
    
    fixer.removeHeterogens()
        
    # Remove all hydrogrens
    hydrogens = [atom for atom in fixer.topology.atoms() if atom.name=='H']
    if len(hydrogens) > 0:
        mod = Modeller(fixer.topology, fixer.positions)
        mod.delete(hydrogens)
        fixer.topology = mod.topology
        fixer.positions = mod.positions

    with open(out_file, 'w') as out_f:
        PDBFile.writeFile(fixer.topology, fixer.positions, out_f, keepIds=True)
        

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="PDBFixer script to add missing atoms")
    parser.add_argument('-i', '--input', help="input file location", required=True)
    parser.add_argument('-o', '--output', help="output file location", required=True)
    args = vars(parser.parse_args())
    fix_file(args['input'], args['output'])