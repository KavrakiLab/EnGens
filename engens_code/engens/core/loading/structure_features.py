"""
Module is responsible for holding all information needed for creating feature readers.
"""
import abc
from abc import ABC
from typing import List, Optional
from fileTypes.interface import SimulationFile
from loading.const import FeatureCharacteristics


class StructureFeatures(ABC):
    """
    Object to contain feature reader which gives structure features to the user
    """

    def __init__(self, ) -> None:
        self._reader: Optional[pyemma.coors.data.FeatureReader] = None

    @classmethod
    @abc.abstractmethod
    def feature_characteristics(cls) -> FeatureCharacteristics:
        """

        :return:
        """
        raise NotImplementedError

    @abc.abstractmethod
    def _routine(self, *args, **kwargs) -> None:
        """
        Generic function for routines needed prior to extracting features.

        :param args:
        :param kwargs:
        :return:
        """
        raise NotImplementedError

    @property
    def reader(self) -> pyemma.coors.data.FeatureReader:
        """
        Feature reader getter

        :return:
        """
        if self._reader is None:
            raise ValueError("Reader has not been set and returned None")
        if not isinstance(self._reader, pyemma.coors.data.FeatureReader):
            raise TypeError("Reader did not have the expected type.")
        return self._reader

    def load(self, structure_source_files: List[SimulationFile], *args, **kwargs) -> None:
        """
        Feature reader setter

        :param structure_source_files: Structure files needed to create feature reader.
        :param args: Arguments needed for feature reader
        :param kwargs: key value arguments needed for feature reader
        :return: Feature reader
        """
        if self._reader is not None:
            raise ValueError(f"Reader appears to already be set. The current reader value is {self._reader}")
        self._routine(*args, **kwargs)
        self._reader = pyemma.coors.data.FeatureReader([x.path for x in structure_source_files], *args, **kwargs)


class StandAloneReader(StructureFeatures):
    """
    No actions are needed prior to feature reader creation. Routine should be empty.
    """

    @classmethod
    def feature_characteristics(cls) -> FeatureCharacteristics:
        return FeatureCharacteristics.StandAlone

    def _routine(self, *args, **kwargs) -> None:
        """
        method is empty purposefully for

        :return:
        """
        return None


class AlignmentBasedReader(StructureFeatures):
    """
    Feature reader that require alignment prior to instantiation of feature reader
    """
    @classmethod
    def feature_characteristics(cls) -> FeatureCharacteristics:
        return FeatureCharacteristics.Aligned

    def _routine(self, *args, **kwargs) -> None:
        pass


class SubstructureBasedReader(StructureFeatures):
    """
    Feature reader that require substructure prior to instantiation of feature reader
    """
    @classmethod
    def feature_characteristics(cls) -> FeatureCharacteristics:
        return FeatureCharacteristics.Aligned

    def _routine(self, *args, **kwargs) -> None:
        pass


def get_feature_reader(feature_characteristic: FeatureCharacteristics) -> StructureFeatures:
    """
    Function loads appropriate feature reader class.

    :param feature_characteristic:
    :return:
    """
    for subclass in StructureFeatures.__subclasses__():
        if subclass.feature_characteristics() == feature_characteristic:
            return subclass()
