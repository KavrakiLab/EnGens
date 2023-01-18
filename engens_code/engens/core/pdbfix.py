#!/miniconda3/bin/python

from openmm.app import PDBFile
from openmm.app import Modeller
from pdbfixer import PDBFixer
import sys
import argparse
import copy

def pdbfix(args):
    input_file = args['input']
    fixer = PDBFixer(filename=input_file)
    output_file = args['output']
    # Find missing residues
    fixer.findMissingResidues()
    chains = list(fixer.topology.chains())
    keys = list(fixer.missingResidues.keys())
    if not args['add_tails']:
        # Do not mess with the residues at the ends of the chains
        for key in keys:
            chain = chains[key[0]]
            if key[1] == 0 or key[1] == len(list(chain.residues())):
                del fixer.missingResidues[key]
            else:
                # do not add wrong long loops
                if len(fixer.missingResidues[key]) > 5:
                    del fixer.missingResidues[key]
    fixer.missingResidues = []
    fixer.findNonstandardResidues()
    if len(fixer.nonstandardResidues) > 0:
        print("Warning: structure has nonstandard residues!")
        print(fixer.nonstandardResidues)
    if args['replace_nonstd']:
        fixer.replaceNonstandardResidues()
    if args["remove_het"] and not args["keep_water"]:
        fixer.removeHeterogens(args["keep_water"])
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    if args['remove_h']:
        # Remove all hydrogrens
        hydrogens = [atom for atom in fixer.topology.atoms() if atom.name=='H']
        if len(hydrogens) > 0:
            mod = Modeller(fixer.topology, fixer.positions)
            mod.delete(hydrogens)
            fixer.topology = mod.topology
            fixer.positions = mod.positions

        with open(output_file, 'w') as out_f:
            PDBFile.writeFile(fixer.topology, fixer.positions, out_f, keepIds=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="PDBFixer script to prepare the receptor for MD.")
    parser.add_argument('-i', '--input', help="input file location", required=True)
    parser.add_argument('-o', '--output', help="output file location", required=True)
    parser.add_argument('--remove_h', help="do you wand pdbfixer to remove original hydrogens?", default=True)
    parser.add_argument('--add_tails', help="do you wand pdbfixer to add missing ends of chains?", default=False)
    parser.add_argument('--replace_nonstd', help="do you wand pdbfixer to replace nonstandard residues?", default=False)
    parser.add_argument('--keep_water', help="do you wand pdbfixer to keep water from the crystal?", default=False)
    parser.add_argument('--remove_het', help="do you wand pdbfixer to remove HET atoms from the crystal?", default=True)

    args = vars(parser.parse_args())
    pdbfix(args)