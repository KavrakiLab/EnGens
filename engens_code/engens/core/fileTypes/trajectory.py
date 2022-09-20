"""
Module is responsible for holding all relevant functionality for trajectory files.
"""
from pathlib import Path
from typing import List

from const import TrajectoryFileTypes
from interface import SimulationFile
from topology import TopologyFile


class TrajectoryFile(SimulationFile):
    """
    Trajectory File class to hold all trajectory files relevant
    """

    def __init__(self, path: Path, topology_file: TopologyFile):
        super().__init__(path)
        self.topology_file: TopologyFile = topology_file

    @classmethod
    def supported_extensions(cls) -> List[str]:
        return TrajectoryFileTypes.get_all()
