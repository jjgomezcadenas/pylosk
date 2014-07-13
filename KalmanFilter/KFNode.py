import sys
from KFBase import KFVector, KFMatrix
from KFMeasurement import KFMeasurement 
from KFZSlice import KFZSlice 
from math import sqrt,log


class KFNode(object):
    """

    A pair (KMeasurement, KZSlice)
    """
    
    def __init__(self,measurement,zslice):
        """
            
        """
        if isinstance(measurement, KFMeasurement)!= True:
            print "The measurement must take an instance of KFMeasurement"
            sys.exit(-1)

        if isinstance(zslice, KFZSlice)!= True:
            print "The zslice must take an instance of KFZSlice"
            sys.exit(-1)

        self.Measurement = measurement 
        self.ZSlice =zslice

    def __str__(self):
        s ="<\nMeasurement :"
        s+=self.Measurement.__str__()+"\nZSlice :"+self.ZSlice.__str__()+"\n>"
        return s

    def __repr__(self):
        return self.__str__() 
    
    
    
if __name__ == '__main__':

    LrXe = 15300 # radiation length of xenon in mm
    Pr = 10      #pressure in bar

    zslice=KFZSlice(LrXe/Pr,101,102,0.1)
    print "KFZSlice =",zslice

    x0=5.
    y0=15.
    hit =KFVector([x0,y0])
    sigma_x=0.1
    sigma_y=0.2
    cov=KFMatrix([[sigma_x**2,0],
                   [0,sigma_y**2]])

    print "hit =",hit

    print "cov =",cov

    measurement = KFMeasurement(hit,cov)

    node =KFNode(measurement,zslice)

    print "node =",node
    
    
    
