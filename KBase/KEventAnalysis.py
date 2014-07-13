"""
An interface to analyse events from file
"""
from abc import ABCMeta, abstractmethod
from KalmanFilter import *

LrXe = 15300. # radiation length of xenon in mm
Pr = 5.      #pressure in bar
sample5 =5.


class KEventAnalysis(object):
    __metaclass__ = ABCMeta

    """
    KEventAnalysis provides an interface to analyse events 
    generated with Irene/Paolina/Toy etc from file.  
    JJ, Spring, 2014

    eventReader --> Implements class KEventReader
    KFitter--->Implements class KalmanFilter
    smearVector --> [sigma_x, sigma_y, sigma_z]
    """

    def __init__(self,eventReader, KFitter, smearVector,
                radiationLength=LrXe,sample=sample5, Pressure=Pr,
                setupType='True'):
    	"""
    	Constructor
    	"""
    	self.eventReader=eventReader
        self.KFitter=KFitter
        self.smearVector=smearVector
        self.radiationLength=radiationLength
        self.sample=sample
        self.Pressure=Pressure
        self.setupType=setupType
        self.betaMax = self.GetKMCParticle()

        self.Hits =[]

        if setupType == 'True':
            self.Hits = self.betaMax.TrueHits()
        else:
            self.Hits = self.betaMax.SmearedHits()


    @abstractmethod
    def GetKMCParticle(self):
        """
        Returns the KMC particle   
        """
        return

    def Hits(self):
        """
        Returns the selected Hits   
        """
        return self.Hits

    def KSystem(self):
        """
        Returns a KFSetup using the hits of the GetKMCParticle 
        setupType = True creates the setup with the true hits
        setupType = Smear creates the setup with the smeared hits 
        """
        
    	Nodes=[]
    	
    	numberOfHits = len(self.Hits)

    	for k in range(0,numberOfHits):
        hitk =self.Hits[k]
        hitkp1 =self.Hits[k+1]
        zki =hitk[2]
        zkf = hitkp1[2]
        edepk = hitk[3]
        xk = hitk[0]
        yk =hitk[1]

        zslice=KFZSlice(radiationLength,zki,zkf,edepk)
        hit =KFVector([xk,yk])
        cov=KFMatrix([[ip.sigma_x**2,0],
                   [0,ip.sigma_y**2]])
                
        measurement = KFMeasurement(hit,cov)
        node =KFNode(measurement,zslice)
        Nodes.append(node)

        setup = KFSetup(Nodes) 
        system =KFSystem(setup)
        return system

    
    def KInitialState(self):
        """
        Returns the initial state using the hits of the GetKMCParticle 
        and the model 
        """

        if self.KFitter.Fitter() == "KFLineal":
            return self.__KFlinealInitState()
        else:
            print "not implemented"
            sys.exit(-1)

    def __KFLinealInitState():
        """
        Initial State for fitter KFLineal 
        """
        hit0 = self.Hits[0]
        x0 = hit0[0]
        y0 = hit0[1]
        ux0,uy0 = self.betaMax.DirectionTangents()
        
        sx0 = self.smearVector[0]
        sy0 = self.smearVector[1]
        sux0 =0.5*ux0 # arbitraryly large relative value
        suy0 =0.5*uxy

        x0x0 = sx0**2
        y0y0 = sy0**2     
        ux0ux0 =sux0**2
        uy0uy0 =suy0**2
        x0ux0 = sx0*sux0
        x0uy0 = sx0*suy0
        y0ux0 = sy0*sux0
        y0uy0 = sy0*suy0
        x0y0 = sx0*sy0
        ux0uy0 = sux0*suy0

        traj0 = KFVector([x0,y0,ux0,uy0])
        cov0 = KFMatrix([
                    [x0x0,x0y0,x0ux0,x0uy0],
                    [x0y0,y0y0,y0ux0,y0uy0],
                    [x0ux0,y0ux0,ux0ux0,ux0uy0],
                    [x0uy0,y0uy0,ux0uy0,uy0uy0]
                    ])


        return (traj0,cov0) 

    
    	