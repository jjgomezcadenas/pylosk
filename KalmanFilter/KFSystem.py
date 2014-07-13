import sys
from KFBase import KFVector, KFMatrix
from KFMeasurement import KFMeasurement 
from KFZSlice import KFZSlice
from KFNode import KFNode 
from KFState import KFState
from KFSetup import KFSetup 
from math import *


class KFSystem(object):
    """

    Represents the system in which the KF operates, including:
    A setup (KFSetup) containing a collection of nodes (measurements, zslices)
    and a container of states (KFState) which describe the states of the trajectory
    in the nodes: 
    
    """
    
    def __init__(self,setup):
        """
            self.Nodes is a list of nodes
            self.States is a list of (dictionary of) states
        """
        if isinstance(setup, KFSetup) == True:
            self.Nodes = setup.Nodes
        elif isinstance(setup, list) == True:
            for node in setup:
                if isinstance(node, KFNode) == False:
                    print "List of nodes must be of type KNode"
                    sys.exit(-1)
            self.Nodes = setup
        else:
            print "the Setup must be a KFSetup or a list of KNode"
            sys.exit(-1)

      
        self.States = []
    

    def __str__(self):
        s ="<System\n Setup: "+self.Nodes.__str__()+"\nStates :"+self.States.__str__()
        s+=">"
        return s
    def __repr__(self):
        return self.__str__() 
    

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

    LrXe = 15300 # radiation length of xenon in mm
    Pr = 10      #pressure in bar

    Nodes=[]
    for i in range(1,6):
        zslice=KFZSlice(LrXe/Pr,100+i,200+i,0.1)
    
        x0=5+i
        y0=15+i
        hit =KFVector([x0,y0])
        sigma_x=0.1
        sigma_y=0.2
        cov=KFMatrix([[sigma_x**2,0],
                   [0,sigma_y**2]])

        measurement = KFMeasurement(hit,cov)

        node =KFNode(measurement,zslice)
        Nodes.append(node)

    setup = KFSetup(Nodes)
    print "Setup =",setup

    traj = StraightLineTrajectory(0,0,1,0.1,0.2)
    Q = StraightLineQ(0,0,1,0.1,0.2)
    
    state = KFState(traj,Q)
    state.Status = "Init"

    print "State=",state

    system =KFSystem(setup)
    system.States.append(state)

    print system

    
    
