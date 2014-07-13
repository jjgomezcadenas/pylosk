"""
Irene Event Reader: Implements the interface KEventReader
"""
from ROOT import *
from KEventReader import KEventReader

gSystem.Load("libirene")

class IEventReader(KEventReader):
    """
    Implements the interface KEventReader for Irene 
    JJ, Spring, 2014
    """

    def __init__(self,pathToFile,gsEvtTree = "EVENT",gsEvtBranch = "EventBranch"):
    	"""
    	Constructor
    	"""

        KEventReader.__init__(self,pathToFile)
    	
        self.fFile = TFile.Open(self.pathToFile)
        self.fEvtTree = self.fFile.Get(gsEvtTree)

        self.ievt =irene.Event()  #create an irene event
        self.fEvtTree.SetBranchAddress(gsEvtBranch, self.ievt)
        self.numberOfBytesRead = 0
        self.totalNumberOfBytesRead =0


    def NumberOfEvents(self):
        """
        Returns the number of events in file  
        """
        return self.fEvtTree.GetEntries()

    def ReadEvent(self,eventNumber):
        """
        Reads event number and returns a handle to the event  
        """
        self.eventNumber = eventNumber
        self.numberOfBytesRead= self.fEvtTree.GetEntry(eventNumber)
        self.totalNumberOfBytesRead+=self.numberOfBytesRead

        return self.ievt  #return handle to irene event

    def NumberOfBytesRead(self):
        return self.numberOfBytesRead

    def TotalNumberOfBytesRead(self):
        return self.totalNumberOfBytesRead

    def EventNumber(self):
        return self.eventNumber

    def CloseFile(self):
        """
        Closes the file  
        """
        self.fFile.Close();
    	