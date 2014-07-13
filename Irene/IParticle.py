"""
Provides a set of wrapper classes to handle Irene objects
JJ, Spring, 2014
"""

import numpy as np
from KMCParticle import KMCParticle



class IParticle(KMCParticle):
    """
    Implements the interface to KMCParticle with an Irene particle
    """

    def __init__(self,ievt,smearVector,sample=5,TMAX=2447,npart=1):
        """
        Grab the main properties of the particle
        TMAX is the Qbb of Xe-136 in keV
        """

        self.TMAX = TMAX
        self.ievt = ievt
        ipart,itrk = self.__SelectEMax(npart)
        self.ihits = itrk.GetHits()
        
        i_vx = ipart.GetInitialVertex()
        i_p = ipart.GetInitialMomentum()
        d_vx = ipart.GetDecayVertex()

        self.name = ipart.Name()
        self.pdg = ipart.GetPDGcode()
        self.mass =ipart.GetMass()
        self.E0 =ipart.Energy()
        self.P0= ipart.Momentum()
        self.T0 =self.E0 - self.mass 
        self.V0 = np.array([i_vx.X(),i_vx.Y(),i_vx.Z()])
        self.Vf = np.array([d_vx.X(),d_vx.Y(),d_vx.Z()])
        self.P30 = np.array([i_p.X(),i_p.Y(),i_p.Z()])
        self.trackLength = ipart.GetTrackLength()
        if(self.P30[2] > 0):
            self.ux = self.P30[0]/self.P30[2]
            self.uy = self.P30[1]/self.P30[2]
        else:
            self.ux = self.P30[0]/self.P0
            self.uy = self.P30[1]/self.P0

        self.THits = self.__FillTrueHits()

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


    def __FillTrueHits(self):
        """
        Fills the true hits in an array
        """

        TrueHits=[]        
        
        for i in range(0,self.ihits.size()):
            ihit = self.ihits.at(i)
            xyzt = ihit.first
            energy = ihit.second 
            hit =np.array([xyzt.X(),xyzt.Y(),xyzt.Z(),energy])
            TrueHits.append(hit)

        return TrueHits

    def __SelectEMax(self,npart):
        """
        Selects the electron that is the (npart)th most energetic
        in the list of particles for the event, where npart = 1 or 2.
        """
    
        eps = 0.01 #tolerance

        itrks= self.ievt.GetTracks()
        n_itrks = itrks.GetEntries()
        
        Emax = -1.; E2max = -1.
        ipmax = 0; ip2max = 0
        itrkmax=0; itrk2max = 0
        for it in range(0,n_itrks):
            itrk = itrks.At(it) 
            ipart= itrk.GetParticle()
            
            if ipart.IsPrimary() == False:
                continue

            T= (ipart.Energy()-ipart.GetMass())*1e+3
            
            # Assign the energies, particles, and tracks correctly.
            if(T > E2max): 
                E2max = T
                ip2max = ipart
                itrk2max = itrk
            if(E2max > Emax):
                Etemp = Emax; Emax = E2max; E2max = Etemp
                iptemp = ipmax; ipmax = ip2max; ip2max = iptemp
                itrktemp = itrkmax; itrkmax = itrk2max; itrk2max = itrktemp
                
#            if abs(T-self.TMAX)< eps:
#                ipmax = ipart
#                itrkmax=itrk
        if(npart == 1):
            return ipmax,itrkmax
        if(npart == 2):
            return ip2max,itrk2max
