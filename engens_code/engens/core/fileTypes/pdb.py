"""
Module is responsible for holding all relevant functionality for PDB files.
"""
from pathlib import Path
from typing import List
from Bio import PDB

from const import PDBFileTypes
from interface import Alignable, SimulationFile, PDBParserTool


class PDBFile(SimulationFile, Alignable):
    """
    Trajectory File class to hold all trajectory files relevant
    """
    def __init__(self, pdb_tool: PDB, path: Path) -> None:
        super(SimulationFile).__init__(path)
        super(PDBFile).__init__(pdb_tool)

    @classmethod
    def supported_extensions(cls) -> List[str]:
        return PDBFileTypes.get_all()

    def align(self, reference_file: SimulationFile, pdb_parser: PDBParserTool ) -> None:

        # first structure as reference structure
        ref_structure = pdb_parser.tool.get_structure("tmp_ref", self.path)
        ref_atoms = []  # only align C- alpha
        # Iterate of all chains in the model in order to find all residues
        for ref_chain in ref_structure[0]:
            # Iterate of all residues in each model in order to find proper atoms
            for ref_res in ref_chain:
                ref_atoms.append(ref_res['CA'])

        x = tqdm.tqdm(self.path, f"Aligning pdb file {self.path}")
        sample_structure = pdb_parser.tool.get_structure("tmp_sample", x)
        sample_model = sample_structure[0]
        sample_atoms = []  # only Calpha
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
        io.save(x[:-4] + "_algn.pdb")
        self.pdb_files = [f[:-4] + "_algn.pdb" for f in self.pdb_files]