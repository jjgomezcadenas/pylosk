import sys
from KFBase import KFVector, KFMatrix
from KFMeasurement import KFMeasurement 
from KFZSlice import KFZSlice
from KFNode import KFNode 
from math import sqrt,log


class KFSetup(object):
    """

    KFSetup is a container of nodes
    
    """
    
    def __init__(self,listOfNodes):
        """
        """
        
        if isinstance(listOfNodes, list)!= True:
            if isinstance(listOfNodes, tuple)!= True:
                print "KFSetup takes a list or tuple of nodes"
            sys.exit(-1)

        for node in listOfNodes:
            if isinstance(node, KFNode)!= True:
                print "nodes must be instance of KFNode"
                sys.exit(-1)


        self.Nodes = listOfNodes

    def __str__(self):
        s ="<\n"
        for node in self.Nodes:
            s+=node.__str__()
        s+="\n>"
        return s
    def __repr__(self):
        return self.__str__() 
    
    
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
    
    
