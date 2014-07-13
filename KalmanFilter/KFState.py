import sys
from KFBase import KFVector, KFMatrix, KFHVector
from math import *

class KFState(KFHVector):
    """

    Represents a State characterized by a trajectory (for example (x,y,ux,uy)) and
    a Covariance matrix. 
    
    """
    
    def __init__(self,traj, cov):
        """
        traj: a KFVector describing the trajectory
        cov a KFMatrix of Covariance 
        
        """

        KFHVector.__init__(self,traj,cov)
        self.Chi2 = 1e+6
    

    def __str__(self):
        return self.__repr__() 

    def __repr__(self):
        s ="<\nChi2 = {0} \n Traj={1} \nCov ={2}\n>".format(
            self.Chi2,self.V.__str__(),self.Cov.__str__())
        return s  

def StraightLineTrajectory(x0,y0,z0,thetax,thetay):
    return KFVector([x0-z0*tan(thetax),
                     y0-z0*tan(thetay),
                     tan(thetax),
                     tan(thetay)
                    ])

def StraightLineQ(x0,y0,z0,thetax,thetay):
    p1 = x0-z0*tan(thetax)
    p2 = y0-z0*tan(thetay)
    p3 = tan(thetax)
    p4 = tan(thetay)
    dz=z0
    dz2=dz*dz
    s2tms =1.

    p3p3 =s2tms*(1+p3**2)*(1+p3**2+p4**2)
    p4p4 = s2tms*(1+p4**2)*(1+p3**2+p4**2)
    p3p4=s2tms*p3*p4*(1+p3**2+p4**2)

    cov=KFMatrix([
        [dz2*p3p3,dz2*p3p4,-dz*p3p3,-dz*p3p4],
        [dz2*p3p4,dz2*p4p4,-dz*p3p4,-dz*p4p4],
        [-dz*p3p3,-dz*p3p4,p3p3,p3p4],
        [-dz*p3p4,-dz*p4p4,p3p4,p4p4]
                   ])

    return cov

    
if __name__ == '__main__':

    traj = StraightLineTrajectory(0,0,1,0.1,0.2)
    Q = StraightLineQ(0,0,1,0.1,0.2)
    
    state = KFState(traj,Q)

    print "State=",state
    
    
    
