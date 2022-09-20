"""
Module is responsible for holding all relevant functionality for topology files.
"""
from typing import List

from interface import SimulationFile
from const import TopologyFileTypes


class TopologyFile(SimulationFile):
    """
    Topology File class to hold all trajectory files relevant

    """

    @classmethod
    def supported_extensions(cls) -> List[str]:
        return TopologyFileTypes.get_all()
