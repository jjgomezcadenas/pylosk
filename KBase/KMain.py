"""
Main driver: a file with event and calls the user-supplied method: 
"""

import sys,os

from abc import ABCMeta, abstractmethod


class KMain(object):
    __metaclass__ = ABCMeta

    """
    KMain is the base class for the driver  
    JJ, Spring, 2014
    """

    def __init__(self,pathToFile,events,generator):
        """
        Constructor: pathToFile is absolute, events in the number
        of events to analyze 
        """

        self.pathToFile=pathToFile
        self.events = events
        self.generator = generator
       
        
    @abstractmethod
    def GetEventReader():
        """
        Returns the handle to the event
        """
        return 
    

