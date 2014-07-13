"""
ToyParticle is a particle whose track was generated with a multiple
scattering toy Monte Carlo
Josh, Spring, 2014
"""
from KMCParticle import KMCParticle

from math import *
from abc import ABCMeta, abstractmethod
import random
import numpy as np

class ToyParticle(KMCParticle):
    __metaclass__ = ABCMeta

    """
    ToyParticle is an electron whose track was generated with a multiple
    scattering toy Monte Carlo
    Josh, Spring, 2014
    """

    def __init__(self,tfile,rev,smearVector,sample=5):
        """
        Initialize the particle
        tfile: the name of the file containing the track
        
        The particle is assumed to be an electron
        Pressure is measured in bar and is used to scale histograms
        """

        self.reverse = rev;   # True if this is a reversed track
        self.name = "electron"
        self.pdg = 11
        self.mass = 0.5109989

        # Generate the track; note this sets T0, E0, T0, P30, ux, uy, 
        #  V0, Vf, trackLength, and TrueHits
        self.THits = self.__ReadTrack(tfile)

        KMCParticle.__init__(self,smearVector,sample)
        

    def Name(self):
        """
        The name of the particle 
        """
        return self.name

    
    def Mass(self):
        """
        The mass of the particle in MeV
        """
        return self.mass

    
    def PDG(self):
        """
        The PDG of the particle 
        """
        return self.pdg

    
    def Energy(self):
        """
        The energy of the particle in MeV at origin
        """
        return self.E0

    
    def KineticEnergy(self):
        """
        The energy of the particle in MeV at origin
        """
        return self.T0

    
    def Momentum(self):
        """
        The momenum of the particle in MeV at origin
        """
        return self.P0

    
    def P3(self):
        """
        The three-momentum of the particle in MeV at origin
        """
        return self.P30

    
    def Vertex(self):
        """
        The production vertex in cm
        """
        return self.V0

    
    def DecayVertex(self):
        """
        The end vertex in cm
        """
        return self.Vf

    def TrackLength(self):
        """
        The track length of this particle in the medium
        """
        return self.trackLength

    
    def DirectionTangents(self):
        """
        ux = Px/Pz and uy = Py/Pz
        """
        return (self.ux,self.uy)

       
    def PropagateToZ(self,z):
        """
        It propagates the trajectory of the particle (assumed a straight line, n MS, no B) to z
        """
        t = (z-self.V0[2])/self.P30[2]
        x = self.V0[0]+self.P30[0]*t
        y = self.V0[1]+self.P30[1]*t
        return x,y

    def TrueHits(self):
        """
        True hits left by this particle in the detector. The particle trajectory could have been
        affected by MS, Eloss and the effect of B, but there is no detector resolution.
        hit = np.array(x,y,z,edep)
        """
        return self.THits
        
    def __ReadTrack(self,tfile):
        """
        Reads in the track from the specified file tfile and constructs
        the corresponding TrueHits array.
        """
        
        TrueHits = []
        
        # Read in the track file.
        # x0 y0 zi zf ux uy uz E deltaE deltaX
        trktbl = np.loadtxt(tfile);
        trk_x0 = trktbl[:,0];
        trk_y0 = trktbl[:,1];
        trk_zi = trktbl[:,2];
        trk_zf = trktbl[:,3];
        trk_ux = trktbl[:,4];
        trk_uy = trktbl[:,5];
        trk_uz = trktbl[:,6];
        trk_E = trktbl[:,7];
        trk_deltaE = trktbl[:,8];
        trk_deltaX = trktbl[:,9];
        
        # Set the initial variables.
        if(self.reverse):
            self.T0 = trk_E[-1];
        else:
            self.T0 = trk_E[0]
        self.E0 = self.T0 + self.mass
        self.P0 = self.E0**2 - self.mass**2
        self.P30 = np.array([0,0,self.P0])
        self.ux = trk_ux[0]
        self.uy = trk_uy[0]

        x0 = trk_x0[0]; y0 = trk_y0[0]; zi0 = trk_zi[0]; zf0 = trk_zf[0];
        xf = trk_x0[-1]; yf = trk_y0[-1]; zif = trk_zi[-1]; zff = trk_zf[-1];
        
        self.V0 = np.array([x0,y0,(zi0+zf0)/2.])
        self.Vf = np.array([xf,yf,(zif+zff)/2.])
        
        # Construct the TrueHits list.
        tlen = 0.
        for x0,y0,zi,zf,ux,uy,uz,deltaE,deltaX in zip(trk_x0,trk_y0,trk_zi,trk_zf,trk_ux,trk_uy,trk_uz,trk_deltaE,trk_deltaX):
            hit = np.array([x0,y0,zi,deltaE])
            tlen += deltaX*sqrt(1 + (uz*ux)**2 + (uz*uy)**2)
            TrueHits.append(hit)
        
        self.trackLength = tlen
        
        return TrueHits;