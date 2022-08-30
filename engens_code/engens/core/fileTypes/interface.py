"""
Module holds all interfaces relevant to file types.
"""
import abc
from abc import ABC
from pathlib import Path
from typing import List, TypeVar

from Bio import PDB

T = TypeVar("T")


class FileTypeNotSupported(Exception):
    """
    Exception class for file types not supported.
    """


class SimulationFile(ABC):
    """
    Abstract class providing common structure for all simulation file types such as trajectory, and PDB files.
    """

    def __init__(self, path: Path):
        """
        To qualify as a simulation file, a path must be provided at the very least.

        :param path: Pure Path of the simulation file
        """
        self.path: Path = path

    @property
    def path(self) -> Path:
        try:
            if self._path is None:
                raise ValueError(f"Paths was not set.")
        except AttributeError:
            raise RuntimeError(f"Path variable was not initialized "
                               f"internally from called class leading to attribute error")
        return self._path

    @path.setter
    def path(self, path: Path) -> None:
        if type(path) != Path:
            raise TypeError(f"Type provided doesn't match what's requested, {type(path)}")
        if not path.exists():
            raise FileNotFoundError(f"Path provided, {path} does not exist.")
        if path.suffix not in SimulationFile.supported_extensions():
            raise FileTypeNotSupported
        self._path = path

    @classmethod
    @abc.abstractmethod
    def supported_extensions(cls) -> List[str]:
        """
        Class method to define supported file types of a given SimulationFile
        :return: List of string representations of allowed file extensions
        """
        raise NotImplementedError


class SimulationTool(ABC):
    """
    A simulation file with actionable tools.
    """
    def __init__(self):
        self._tool: T = None

    @abc.abstractmethod
    @property
    def tool(self) -> T:
        raise NotImplementedError


class PDBParserTool(SimulationTool):
    """
    PDB parser tool
    """

    @property
    def tool(self) -> PDB.PDBParser:
        if self._tool is None:
            self._tool = PDB.PDBParser(QUIET=True)
        return self._tool


class Alignable(ABC):

    def __init__(self) -> None:
        self.aligned: bool = False

    @abc.abstractmethod
    def align(self, simulation_file: SimulationFile) -> None:
        """
        Performs alignment on self. I.E. trajectory file alignment, PDB alignment.

        :param simulation_file: Simulation file needed for alignment such as a reference structure.
        :return: None
        """
        raise NotImplementedError
