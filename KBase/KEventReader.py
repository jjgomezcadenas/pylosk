"""
An interface to read events from file
"""
from abc import ABCMeta, abstractmethod

class KEventReader(object):
    __metaclass__ = ABCMeta

    """
    KEventReader provides an interface to read events 
    generated with Irene/Paolina/Toy etc from file.  
    JJ, Spring, 2014
    """

    def __init__(self,pathToFile):
    	"""
    	Constructor
    	"""
    	self.pathToFile=pathToFile

    @abstractmethod
    def NumberOfEvents(self):
        """
        Returns the number of events in file  
        """
        return

    @abstractmethod
    def ReadEvent(self,eventNumber):
        """
        Reads event number and returns a handle to the event  
        """
        return

    @abstractmethod
    def EventNumber():
        """
        Returns the event number 
        """
        return

    @abstractmethod
    def CloseFile(self):
        """
        Closes the file  
        """
        return
    	