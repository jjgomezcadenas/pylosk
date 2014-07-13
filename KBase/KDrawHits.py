"""
Base class to draw hits
"""

from abc import ABCMeta, abstractmethod
import numpy as np

class KDrawHits(object):
    __metaclass__ = ABCMeta
    """
    Define an interface to draw hits & measurements together.
    """

    def __init__(self,binLimits,measurements,hits):
        """
        binLimits =([xl,yl,zl],[xu,yu,zu])
        measurements = [x,y,z,e] or numpy array
        """

        self.lbin = binLimits[0]
        self.ubin = binLimits[1]
        self.Meas = measurements
        self.Hits = hits


    @abstractmethod 
    def DrawMeasurements(self, draw='2D'):
        """
        Draws Measurements: the field draw selects 2D or 3D.
        """

    @abstractmethod 
    def DrawHits(self, draw='2D'):
        """
        Draws hits: the field draw selects 2D or 3D.
        """

    @abstractmethod 
    def DrawAll(self, draw='2D'):
        """
        Draws measurements & hits: the field draw selects 2D or 3D.
        """
    