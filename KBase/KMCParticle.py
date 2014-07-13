"""
KMCParticle represents a MC particle. It is a base class which can be implemented
with concrete particles (such as IParticle to describe a MC Irene particle)
JJ, Spring, 2014
"""

from abc import ABCMeta, abstractmethod
import random
import numpy as np

class KMCParticle(object):
    __metaclass__ = ABCMeta

    """
    KMCParticle represents a MC particle. It is a base class which can be implemented
    with concrete particles (such as IParticle to describe a MC Irene particle)
    JJ, Spring, 2014
    """

    def __init__(self,smearVector,sample=5):

        self.smearVector =smearVector
        self.sample = sample

        self.SmearHits=self.__FillSmearHits()

    @abstractmethod
    def Name(self):
        """
        The name of the particle 
        """
        return

    @abstractmethod
    def Mass(self):
        """
        The mass of the particle in MeV
        """
        return

    @abstractmethod
    def PDG(self):
        """
        The PDG of the particle 
        """
        return

    @abstractmethod
    def Energy(self):
        """
        The energy of the particle in MeV at origin
        """
        return 

    @abstractmethod
    def KineticEnergy(self):
        """
        The energy of the particle in MeV at origin
        """
        return

    @abstractmethod
    def Momentum(self):
        """
        The momenum of the particle in MeV at origin
        """
        return

    @abstractmethod
    def P3(self):
        """
        The three-momentum of the particle in MeV at origin
        """
        return 

    @abstractmethod
    def Vertex(self):
        """
        The production vertex in cm
        """
        return 

    @abstractmethod
    def DecayVertex(self):
        """
        The end vertex in cm
        """
        return 

    @abstractmethod
    def TrackLength(self):
        """
        The track length of this particle in the medium
        """
        return 

    @abstractmethod
    def DirectionTangents(self):
        """
        ux = Px/Pz and uy = Py/Pz
        """
        return 

    @abstractmethod    
    def PropagateToZ(self,z):
        """
        It propagates the trajectory of the particle (assumed a straight line, n MS, no B) to z
        """
        return 

    @abstractmethod 
    def TrueHits(self):
        """
        True hits left by this particle in the detector. The particle trajectory could have been
        affected by MS, Eloss and the effect of B, but there is no detector resolution.
        hit = np.array(x,y,z,edep)
        """
        return 

    def SmearedHits(self,reverse=False):
        """
        Smear true hits by errors sigma_x,sigma_y,sigma_z.
        Sample the true hits every sample
        """

        # Reverse the list of smeared hits.
        if(reverse):
            rlist = self.SmearHits;
            rlist.reverse();
            return rlist;
        else:
            return self.SmearHits;

    # @abstractmethod 
    # def DrawHits(self, draw='2D',view='True'):
    #     """
    #     Draws hits: the field draw selects 2D or 3D.
    #     The field view can be equal to True for True hits, Smear for Smeared hits or both
    #     """
    #     return

    def __FillSmearHits(self):
        """
        Fills the smeared hits in an array
        """

        SHits=self.__SampleHits()
        SmearHits=[]

        for hit in SHits:
            x = random.gauss(hit[0],self.smearVector[0])
            y = random.gauss(hit[1],self.smearVector[1])
            z = random.gauss(hit[2],self.smearVector[2])
            smhit = np.array([x,y,z,hit[3]])
            SmearHits.append(smhit)

        return SmearHits


    def __SampleHits(self):
        """
        Sample the true hits according to sample
        """

        SHits=[]       
        nhit = 0
        edep=0
        Hits =self.TrueHits()
        n_trueHits = len(Hits)
        
        for i in range(0,n_trueHits):
            hit = Hits[i]
            
            if i==0:
                SHits.append(hit)

            elif nhit==self.sample:
                edep+=hit[3]
                nhit =np.array([hit[0],hit[1],hit[2],edep])
                SHits.append(nhit)
                nhit=0
                edep = 0
            else:
                nhit+=1
                edep+=hit[3]

        return SHits        

    def PrintHits(self,sample,hitType='True'):
        """
        Print True or Smeared Hits (sample them according to sample)
        """

        print self.ViewHits(sample,hitType)

    def ViewHits(self,sample=10,hitType='True'):
        """
        View True or Smeared Hits (sample them according to sample)
        """

        Hits = self.SmearHits
        if hitType=='True':
            Hits = self.TrueHits()

        s="\n"
        for i in range(0,len(Hits)-1,sample):
            hit =Hits[i]
            s+= """
            hit: (x,y,z) =({0},{1},{2}) cm edep ={3} MeV\n
            """.format(hit[0],hit[1],hit[2],hit[3])
        return s


    def __str__(self):
        ux,uy  = self.DirectionTangents()

        s ="""<Name= {0}; PDG ={1}; Mass ={2} MeV;\n
            Energy = {3} MeV; Momentum = {4} MeV; Kinetic Energy = {5} MeV; \n
            Production Vertex ={6} cm; Decay Vertex ={7} cm; Track length ={8} cm; \n
            Initial momentum = {9} MeV;
            ux = tag(theta_x) = tx = {10}
            uy = tag(theta_y) = ty = {11}>
            """.format(
                self.Name(),self.PDG(),self.Mass(),
                self.Energy(),self.Momentum(),self.KineticEnergy(),
                self.Vertex(),self.DecayVertex(),self.TrackLength(),self.P3(),ux,uy
                )
        return s

    def __repr__(self):
        return self.__str__()

