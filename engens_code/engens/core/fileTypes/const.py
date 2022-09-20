"""
Holds all constants for file types.
"""
from enum import Enum
from typing import List


class SupportedFileTypes(Enum):

    @classmethod
    def get_all(cls) -> List:
        file_types = []
        file_types.extend([(x, x.value) for x in cls])
        for file_type in cls.__subclasses__():
            file_types.extend([(x, x.value) for x in file_type])
        return file_types


class TrajectoryFileTypes(SupportedFileTypes):
    XTC = "xtc"
    DCD = "dcd"
    TRR = "trr"
    BINPOS = "binpos"
    NETCDF = "netcdf"
    ARC = "arc"
    HDF5 = "hdf5"
    LAMMPSTRJ = "lammpstrj"


class PDBFileTypes(SupportedFileTypes):
    PDB = "pdb"


class TopologyFileTypes(SupportedFileTypes):
    TPR = "tpr"
    PDB = "pdb"
