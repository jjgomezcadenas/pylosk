
import sys
from KFBase import KFVector, KFMatrix
from math import sqrt,log

class KFZSlice(object):
    """

    A Setup can be sliced along the z axis.
    Every KFZSlice is defined by is thickness, dz = zf-zi,
    and the material that fills it (X0). When an electron
    crosses a ZSlice it will be deflected by a theta_ms
    angle and will loose some energy (Edep)
    """
    
    def __init__(self,Lr,zi,zf,edep):
        """
           Lr is the radiation length of the setup in mm
           zi, zf in mm
           edep in MeV
        """
        

        self.Zi = zi
        self.Zf = zf
        self.dz = zf-zi
        self.Edep = edep
        self.Lr = Lr

    def __str__(self):
        s ="<LR : %7.1f mm, zi : %7.1f mm zf : %7.1f mm edep : %7.1e MeV"%(
            self.Lr, self.Zi,self.Zf,self.Edep)
        s+=">"
        return s
    def __repr__(self):
        return self.__str__() 
    
    
if __name__ == '__main__':

    LrXe = 15300 # radiation length of xenon in mm
    Pr = 10      #pressure in bar


    s=KFZSlice(LrXe/Pr,101,102,0.1)
    print "KZSlice =",s
   
    
    
    
