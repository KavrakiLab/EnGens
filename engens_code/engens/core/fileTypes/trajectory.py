"""
Module is responsible for holding all relevant functionality for trajectory files.
"""
from typing import List

from const import TrajectoryFileTypes
from interface import SimulationFile, Alignable


class TrajectoryFile(SimulationFile, Alignable):
    """
    Trajectory File class to hold all trajectory files relevant
    """

    @classmethod
    def supported_extensions(cls) -> List[str]:
        return TrajectoryFileTypes.get_all()

    def align(self, simulation_file: SimulationFile) -> None:
        pass
