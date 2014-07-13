"""
Base class to draw a vector
"""

from abc import ABCMeta, abstractmethod
import numpy as np

class KDrawVector(object):
    __metaclass__ = ABCMeta
    """
    Define an interface to draw hits & measurements together.
    """

    def __init__(self,vector):
        """
        vector = [x,y,z,e] or numpy array
        """

        self.vector = vector

    @abstractmethod 
    def Draw(self):
        """
        Draws Vector.
        """
