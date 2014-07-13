"""
MatPlotLib Implementation of Base class DrawVector
"""

from KDrawVector import KDrawVector
import matplotlib.pyplot as plt
import os
import sys
import numpy as np

from mpl_toolkits.mplot3d import Axes3D

class MPLDrawVector(KDrawVector):
    """
    Define an interface to draw hits
    """

    def __init__(self,vector):
        """
        vector =[] or np
        """

        KDrawVector.__init__(self,vector)
        
        

    def Draw(self):
        """
        Draws Vector
        """
        
        plt.show(plt.plot(self.vector))
            
        