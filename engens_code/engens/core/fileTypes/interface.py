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


class Loadable(ABC):

    @abc.abstractmethod
    def load(self, *args, **kwargs) -> None:
        """
        Loads file into desired object

        :return: None
        """
        raise NotImplementedError


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
    def path(self, pure_path: Path) -> None:
        if type(pure_path) != Path:
            raise TypeError(f"Type provided doesn't match what's requested, {type(pure_path)}")
        if not pure_path.exists():
            raise FileNotFoundError(f"Path provided, {pure_path} does not exist.")
        if pure_path.suffix not in SimulationFile.supported_extensions():
            raise FileTypeNotSupported
        self._path = pure_path

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
        self._aligned: bool = False

    @property
    def aligned(self) -> bool:
        if self._aligned is None:
            raise ValueError("Alignment value not set.")
        return self._aligned

    @aligned.setter
    def aligned(self, alignment: bool):
        if type(alignment) != bool:
            raise TypeError(f"Did not pass in type, boolean, instead {type(alignment)}.")
        self._aligned = alignment

    @abc.abstractmethod
    def align(self, simulation_file: SimulationFile, tool: SimulationTool, **kwargs) -> None:
        """
        Performs alignment on self. I.E. trajectory file alignment, PDB alignment.

        :param simulation_file: Simulation file needed for alignment such as a reference structure.
        :param tool: tool used to apply alignment. Should be of type simulation tool
        :return: None
        """
        raise NotImplementedError
