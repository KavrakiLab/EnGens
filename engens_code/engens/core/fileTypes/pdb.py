"""
Module is responsible for holding all relevant functionality for PDB files.
"""
from typing import List

from const import PDBFileTypes
from interface import SimulationFile


class PDBFile(SimulationFile):
    """
    Trajectory File class to hold all trajectory files relevant
    """

    @classmethod
    def supported_extensions(cls) -> List[str]:
        return PDBFileTypes.get_all()
